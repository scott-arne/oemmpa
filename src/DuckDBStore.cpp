#include "oemmpa/DuckDBStore.h"

#include "oemmpa/Error.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/PairScoring.h"

#include <duckdb.hpp>

#include <algorithm>
#include <cstdint>
#include <exception>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace OEMMPA {
namespace {

std::string normalize_database_path(const std::string& database_path) {
    return database_path.empty() ? ":memory:" : database_path;
}

const std::vector<std::string>& base_table_names() {
    static const std::vector<std::string> tables = {
        "compound",
        "compound_property",
        "constant_smiles",
        "dataset",
        "environment_fingerprint",
        "pair",
        "property_name",
        "rule",
        "rule_environment",
        "rule_smiles",
    };
    return tables;
}

std::string resolve_table_alias(const std::string& table_name) {
    static const std::unordered_map<std::string, std::string> aliases = {
        {"constants", "constant_smiles"},
        {"molecule_properties", "compound_property"},
        {"molecules", "compound"},
        {"pairs", "pair"},
        {"transforms", "rule"},
    };
    const auto alias_iter = aliases.find(table_name);
    if (alias_iter == aliases.end()) {
        return table_name;
    }
    return alias_iter->second;
}

bool is_base_table_name(const std::string& table_name) {
    const std::string resolved_table_name = resolve_table_alias(table_name);
    const std::vector<std::string>& tables = base_table_names();
    return std::find(tables.begin(), tables.end(), resolved_table_name) != tables.end();
}

duckdb::Value string_or_null(const std::string& value) {
    if (value.empty()) {
        return duckdb::Value(nullptr);
    }
    return duckdb::Value(value);
}

std::string trim_copy(const std::string& text) {
    const std::string whitespace = " \t\r\n";
    const std::string::size_type first = text.find_first_not_of(whitespace);
    if (first == std::string::npos) {
        return "";
    }
    const std::string::size_type last = text.find_last_not_of(whitespace);
    return text.substr(first, last - first + 1);
}

std::string make_generated_external_id(std::uint64_t internal_id) {
    return "molecule_" + std::to_string(internal_id);
}

std::vector<std::string> parse_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    bool in_quotes = false;

    for (std::size_t index = 0; index < line.size(); ++index) {
        const char ch = line[index];
        if (ch == '"') {
            if (in_quotes && index + 1 < line.size() && line[index + 1] == '"') {
                field.push_back('"');
                ++index;
            } else {
                in_quotes = !in_quotes;
            }
        } else if (ch == ',' && !in_quotes) {
            fields.push_back(trim_copy(field));
            field.clear();
        } else {
            field.push_back(ch);
        }
    }

    if (in_quotes) {
        throw StorageError("unterminated quoted CSV field");
    }

    fields.push_back(trim_copy(field));
    return fields;
}

void require_success(const duckdb::QueryResult& result, const std::string& sql) {
    if (result.HasError()) {
        throw StorageError("DuckDB query failed: " + result.GetError() + " SQL: " + sql);
    }
}

std::unique_ptr<duckdb::QueryResult> execute_prepared(
    const std::unique_ptr<duckdb::Connection>& connection,
    const std::string& sql,
    duckdb::vector<duckdb::Value> values
) {
    if (!connection) {
        throw StorageError("DuckDB connection is not open");
    }

    std::unique_ptr<duckdb::PreparedStatement> statement = connection->Prepare(sql);
    if (!statement) {
        throw StorageError("DuckDB prepare returned no statement: " + sql);
    }
    if (statement->HasError()) {
        throw StorageError("DuckDB prepare failed: " + statement->GetError() + " SQL: " + sql);
    }

    std::unique_ptr<duckdb::QueryResult> result = statement->Execute(values);
    if (!result) {
        throw StorageError("DuckDB query returned no result: " + sql);
    }
    require_success(*result, sql);
    return result;
}

std::uint64_t get_next_id(
    const std::unique_ptr<duckdb::Connection>& connection,
    const std::string& table_name,
    const std::string& id_column
) {
    const std::string sql = "select coalesce(max(" + id_column + "), 0) + 1 from " + table_name;
    std::unique_ptr<duckdb::QueryResult> result = connection->Query(sql);
    if (!result) {
        throw StorageError("DuckDB query returned no result: " + sql);
    }
    require_success(*result, sql);

    for (const auto& row : *result) {
        return row.GetValue<std::uint64_t>(0);
    }

    throw StorageError("DuckDB did not return a next ID for " + table_name);
}

std::uint64_t find_named_row_id(
    const std::unique_ptr<duckdb::Connection>& connection,
    const std::string& table_name,
    const std::string& id_column,
    const std::string& value_column,
    const std::string& value
) {
    const std::string sql =
        "select " + id_column + " from " + table_name + " where " + value_column + " = ?";
    duckdb::vector<duckdb::Value> values = {duckdb::Value(value)};
    std::unique_ptr<duckdb::QueryResult> result =
        execute_prepared(connection, sql, std::move(values));

    for (const auto& row : *result) {
        return row.GetValue<std::uint64_t>(0);
    }

    return 0;
}

bool has_external_molecule_id(
    const std::unique_ptr<duckdb::Connection>& connection,
    const std::string& external_id
) {
    const std::string sql = "select count(*) from compound where public_id = ?";
    duckdb::vector<duckdb::Value> values = {duckdb::Value(external_id)};
    std::unique_ptr<duckdb::QueryResult> result =
        execute_prepared(connection, sql, std::move(values));

    for (const auto& row : *result) {
        return row.GetValue<std::int64_t>(0) > 0;
    }
    return false;
}

std::uint64_t find_molecule_internal_id_by_external_id(
    const std::unique_ptr<duckdb::Connection>& connection,
    const std::string& external_id
) {
    const std::string sql = "select id from compound where public_id = ?";
    duckdb::vector<duckdb::Value> values = {duckdb::Value(external_id)};
    std::unique_ptr<duckdb::QueryResult> result =
        execute_prepared(connection, sql, std::move(values));

    for (const auto& row : *result) {
        return row.GetValue<std::uint64_t>(0);
    }
    return 0;
}

std::unordered_map<std::string, std::size_t> build_header_index(
    const std::vector<std::string>& header
) {
    if (header.empty()) {
        throw StorageError("CSV file must contain a header row");
    }

    std::unordered_map<std::string, std::size_t> index_by_name;
    for (std::size_t index = 0; index < header.size(); ++index) {
        const std::string& name = header[index];
        if (name.empty()) {
            throw StorageError("CSV header names must not be blank");
        }
        if (!index_by_name.emplace(name, index).second) {
            throw StorageError("duplicate CSV header: " + name);
        }
    }
    return index_by_name;
}

std::string resolve_id_column(
    const std::unordered_map<std::string, std::size_t>& header_index,
    const std::string& requested_id_column
) {
    if (requested_id_column.empty()) {
        throw StorageError("CSV id column must not be blank");
    }

    if (header_index.find(requested_id_column) != header_index.end()) {
        return requested_id_column;
    }

    if (requested_id_column == "id") {
        for (const std::string& candidate : {"ID", "Name", "name"}) {
            if (header_index.find(candidate) != header_index.end()) {
                return candidate;
            }
        }
    }

    throw StorageError("CSV id column not found: " + requested_id_column);
}

std::vector<std::string> resolve_property_columns(
    const std::vector<std::string>& header,
    const std::unordered_map<std::string, std::size_t>& header_index,
    const std::string& id_column,
    const std::vector<std::string>& requested_property_columns
) {
    std::vector<std::string> property_columns;
    if (requested_property_columns.empty()) {
        for (const std::string& name : header) {
            if (name != id_column) {
                property_columns.push_back(name);
            }
        }
    } else {
        property_columns = requested_property_columns;
    }

    if (property_columns.empty()) {
        throw StorageError("CSV file must contain at least one property column");
    }

    std::set<std::string> seen;
    for (const std::string& name : property_columns) {
        if (name.empty()) {
            throw StorageError("CSV property column names must not be blank");
        }
        if (name == id_column) {
            throw StorageError("CSV id column cannot also be a property column: " + name);
        }
        if (header_index.find(name) == header_index.end()) {
            throw StorageError("CSV property column not found: " + name);
        }
        if (!seen.insert(name).second) {
            throw StorageError("duplicate CSV property column: " + name);
        }
    }
    return property_columns;
}

double parse_property_value(const std::string& property_name, const std::string& text) {
    std::size_t consumed = 0;
    try {
        const double value = std::stod(text, &consumed);
        if (consumed != text.size()) {
            throw std::invalid_argument("trailing characters");
        }
        return value;
    } catch (const std::exception&) {
        throw StorageError(
            "property " + property_name + " value " + text +
            " cannot be converted to a number"
        );
    }
}

std::uint64_t get_or_create_named_row_id(
    const std::unique_ptr<duckdb::Connection>& connection,
    const std::string& table_name,
    const std::string& id_column,
    const std::string& value_column,
    const std::string& value
) {
    if (value.empty()) {
        throw StorageError("cannot store empty normalized value in " + table_name);
    }

    const std::uint64_t existing_id =
        find_named_row_id(connection, table_name, id_column, value_column, value);
    if (existing_id != 0) {
        return existing_id;
    }

    const std::uint64_t next_id = get_next_id(connection, table_name, id_column);
    const std::string sql =
        "insert into " + table_name + " (" + id_column + ", " + value_column + ") values (?, ?)";
    duckdb::vector<duckdb::Value> values = {
        duckdb::Value::UBIGINT(next_id),
        duckdb::Value(value),
    };
    execute_prepared(connection, sql, std::move(values));
    return next_id;
}

std::uint64_t find_rule_id(
    const std::unique_ptr<duckdb::Connection>& connection,
    std::uint64_t from_smiles_id,
    std::uint64_t to_smiles_id
) {
    const std::string sql =
        "select id from rule where from_smiles_id = ? and to_smiles_id = ?";
    duckdb::vector<duckdb::Value> values = {
        duckdb::Value::UBIGINT(from_smiles_id),
        duckdb::Value::UBIGINT(to_smiles_id),
    };
    std::unique_ptr<duckdb::QueryResult> result =
        execute_prepared(connection, sql, std::move(values));

    for (const auto& row : *result) {
        return row.GetValue<std::uint64_t>(0);
    }
    return 0;
}

std::uint64_t get_or_create_rule_id(
    const std::unique_ptr<duckdb::Connection>& connection,
    std::uint64_t from_smiles_id,
    std::uint64_t to_smiles_id
) {
    const std::uint64_t existing_id = find_rule_id(connection, from_smiles_id, to_smiles_id);
    if (existing_id != 0) {
        return existing_id;
    }

    const std::uint64_t next_id = get_next_id(connection, "rule", "id");
    const std::string sql =
        "insert into rule (id, from_smiles_id, to_smiles_id) values (?, ?, ?)";
    duckdb::vector<duckdb::Value> values = {
        duckdb::Value::UBIGINT(next_id),
        duckdb::Value::UBIGINT(from_smiles_id),
        duckdb::Value::UBIGINT(to_smiles_id),
    };
    execute_prepared(connection, sql, std::move(values));
    return next_id;
}

std::uint64_t find_environment_fingerprint_id(
    const std::unique_ptr<duckdb::Connection>& connection,
    const std::string& smarts,
    const std::string& pseudosmiles,
    const std::string& parent_smarts
) {
    const std::string sql =
        "select id from environment_fingerprint "
        "where smarts = ? and pseudosmiles = ? and parent_smarts = ?";
    duckdb::vector<duckdb::Value> values = {
        duckdb::Value(smarts),
        duckdb::Value(pseudosmiles),
        duckdb::Value(parent_smarts),
    };
    std::unique_ptr<duckdb::QueryResult> result =
        execute_prepared(connection, sql, std::move(values));

    for (const auto& row : *result) {
        return row.GetValue<std::uint64_t>(0);
    }
    return 0;
}

std::uint64_t get_or_create_environment_fingerprint_id(
    const std::unique_ptr<duckdb::Connection>& connection,
    const std::string& smarts,
    const std::string& pseudosmiles,
    const std::string& parent_smarts
) {
    const std::uint64_t existing_id = find_environment_fingerprint_id(
        connection,
        smarts,
        pseudosmiles,
        parent_smarts
    );
    if (existing_id != 0) {
        return existing_id;
    }

    const std::uint64_t next_id = get_next_id(connection, "environment_fingerprint", "id");
    const std::string sql =
        "insert into environment_fingerprint "
        "(id, smarts, pseudosmiles, parent_smarts) values (?, ?, ?, ?)";
    duckdb::vector<duckdb::Value> values = {
        duckdb::Value::UBIGINT(next_id),
        duckdb::Value(smarts),
        duckdb::Value(pseudosmiles),
        duckdb::Value(parent_smarts),
    };
    execute_prepared(connection, sql, std::move(values));
    return next_id;
}

std::uint64_t find_rule_environment_id(
    const std::unique_ptr<duckdb::Connection>& connection,
    std::uint64_t rule_id,
    std::uint64_t environment_fingerprint_id,
    int radius
) {
    const std::string sql =
        "select id from rule_environment "
        "where rule_id = ? and environment_fingerprint_id = ? and radius = ?";
    duckdb::vector<duckdb::Value> values = {
        duckdb::Value::UBIGINT(rule_id),
        duckdb::Value::UBIGINT(environment_fingerprint_id),
        duckdb::Value::INTEGER(radius),
    };
    std::unique_ptr<duckdb::QueryResult> result =
        execute_prepared(connection, sql, std::move(values));

    for (const auto& row : *result) {
        return row.GetValue<std::uint64_t>(0);
    }
    return 0;
}

std::uint64_t get_or_create_rule_environment_id(
    const std::unique_ptr<duckdb::Connection>& connection,
    std::uint64_t rule_id,
    std::uint64_t environment_fingerprint_id,
    int radius
) {
    const std::uint64_t existing_id = find_rule_environment_id(
        connection,
        rule_id,
        environment_fingerprint_id,
        radius
    );
    if (existing_id != 0) {
        return existing_id;
    }

    const std::uint64_t next_id = get_next_id(connection, "rule_environment", "id");
    const std::string sql =
        "insert into rule_environment "
        "(id, rule_id, environment_fingerprint_id, radius, num_pairs) "
        "values (?, ?, ?, ?, ?)";
    duckdb::vector<duckdb::Value> values = {
        duckdb::Value::UBIGINT(next_id),
        duckdb::Value::UBIGINT(rule_id),
        duckdb::Value::UBIGINT(environment_fingerprint_id),
        duckdb::Value::INTEGER(radius),
        duckdb::Value::UINTEGER(0),
    };
    execute_prepared(connection, sql, std::move(values));
    return next_id;
}

void increment_rule_environment_pair_count(
    const std::unique_ptr<duckdb::Connection>& connection,
    std::uint64_t rule_environment_id
) {
    const std::string sql =
        "update rule_environment set num_pairs = num_pairs + 1 where id = ?";
    duckdb::vector<duckdb::Value> values = {
        duckdb::Value::UBIGINT(rule_environment_id),
    };
    execute_prepared(connection, sql, std::move(values));
}

bool compare_pairs(const MatchedPair& lhs, const MatchedPair& rhs) {
    return std::make_tuple(
        lhs.GetConstantSmiles(),
        lhs.GetSourceMoleculeId(),
        lhs.GetTargetMoleculeId(),
        lhs.GetTransformSmiles(),
        lhs.GetSourceVariableSmiles(),
        lhs.GetTargetVariableSmiles(),
        lhs.GetCutCount(),
        lhs.GetHeavyAtomDelta(),
        lhs.GetHeavyBondDelta()
    ) < std::make_tuple(
        rhs.GetConstantSmiles(),
        rhs.GetSourceMoleculeId(),
        rhs.GetTargetMoleculeId(),
        rhs.GetTransformSmiles(),
        rhs.GetSourceVariableSmiles(),
        rhs.GetTargetVariableSmiles(),
        rhs.GetCutCount(),
        rhs.GetHeavyAtomDelta(),
        rhs.GetHeavyBondDelta()
    );
}

using PairScoringKey = std::tuple<std::string, unsigned int, unsigned int>;

std::vector<MatchedPair> apply_pair_scoring(
    const std::vector<MatchedPair>& pairs,
    const ScoringOptions& scoring_options
) {
    std::map<PairScoringKey, std::vector<MatchedPair>> candidates_by_group;
    for (const MatchedPair& pair : pairs) {
        candidates_by_group[{
            pair.GetConstantSmiles(),
            pair.GetSourceMoleculeId(),
            pair.GetTargetMoleculeId()
        }].push_back(pair);
    }

    std::vector<MatchedPair> selected_pairs;
    for (const auto& entry : candidates_by_group) {
        std::vector<MatchedPair> selected =
            PairScoring::Select(entry.second, scoring_options);
        selected_pairs.insert(selected_pairs.end(), selected.begin(), selected.end());
    }

    std::sort(selected_pairs.begin(), selected_pairs.end(), compare_pairs);
    return selected_pairs;
}

std::string build_pair_query(const QueryOptions& options) {
    std::ostringstream sql;
    sql
        << "select "
        << "p.compound1_id, p.compound2_id, "
        << "coalesce(source_molecule.public_id, ''), "
        << "coalesce(target_molecule.public_id, ''), "
        << "source_molecule.clean_smiles, target_molecule.clean_smiles, "
        << "c.smiles, "
        << "source_variable.smiles, target_variable.smiles, "
        << "p.cut_count, p.heavy_atom_delta, p.heavy_bond_delta "
        << "from pair p "
        << "join compound source_molecule on source_molecule.id = p.compound1_id "
        << "join compound target_molecule on target_molecule.id = p.compound2_id "
        << "join constant_smiles c on c.id = p.constant_id "
        << "join rule_environment rule_env on rule_env.id = p.rule_environment_id "
        << "join rule r on r.id = rule_env.rule_id "
        << "join rule_smiles source_variable on source_variable.id = r.from_smiles_id "
        << "join rule_smiles target_variable on target_variable.id = r.to_smiles_id "
        << "where true ";

    if (!options.GetSymmetric()) {
        sql << "and p.compound1_id < p.compound2_id ";
    }
    if (options.GetMaxHeavyAtomChange() >= 0) {
        sql << "and abs(p.heavy_atom_delta) <= " << options.GetMaxHeavyAtomChange() << " ";
    }
    if (options.GetMaxRelativeHeavyAtomChange() >= 0.0) {
        sql
            << "and case "
            << "when source_molecule.clean_num_heavies = 0 "
            << "then abs(p.heavy_atom_delta) = 0 "
            << "else cast(abs(p.heavy_atom_delta) as double) / "
            << "cast(source_molecule.clean_num_heavies as double) <= "
            << options.GetMaxRelativeHeavyAtomChange() << " "
            << "end ";
    }

    sql << "order by p.id";
    return sql.str();
}

}  // namespace

DuckDBStore::DuckDBStore()
    : DuckDBStore(":memory:") {}

DuckDBStore::DuckDBStore(const std::string& database_path)
    : database_path_(normalize_database_path(database_path)) {
    try {
        database_ = std::make_unique<duckdb::DuckDB>(database_path_);
        connection_ = std::make_unique<duckdb::Connection>(*database_);
    } catch (const std::exception& exc) {
        throw StorageError("failed to open DuckDB database: " + std::string(exc.what()));
    }
}

DuckDBStore::~DuckDBStore() = default;

void DuckDBStore::InitializeSchema() {
    Execute("begin transaction");
    try {
        Execute(
            "create table if not exists dataset ("
            "id ubigint primary key,"
            "oemmpa_schema_version uinteger not null,"
            "title varchar not null,"
            "creation_date timestamp default current_timestamp,"
            "fragment_options varchar not null,"
            "index_options varchar not null,"
            "is_symmetric boolean not null,"
            "num_compounds ubigint,"
            "num_rules ubigint,"
            "num_pairs ubigint,"
            "num_rule_environments ubigint"
            ")"
        );
        Execute(
            "create table if not exists compound ("
            "id ubigint primary key,"
            "public_id varchar unique,"
            "input_smiles varchar not null,"
            "clean_smiles varchar not null,"
            "clean_num_heavies uinteger not null"
            ")"
        );
        Execute(
            "create table if not exists property_name ("
            "id ubigint primary key,"
            "name varchar not null unique"
            ")"
        );
        Execute(
            "create table if not exists compound_property ("
            "id ubigint primary key,"
            "compound_id ubigint not null references compound(id),"
            "property_name_id ubigint not null references property_name(id),"
            "value double not null,"
            "unique (compound_id, property_name_id)"
            ")"
        );
        Execute(
            "create table if not exists rule_smiles ("
            "id ubigint primary key,"
            "smiles varchar not null unique,"
            "num_heavies uinteger"
            ")"
        );
        Execute(
            "create table if not exists rule ("
            "id ubigint primary key,"
            "from_smiles_id ubigint not null references rule_smiles(id),"
            "to_smiles_id ubigint not null references rule_smiles(id),"
            "unique (from_smiles_id, to_smiles_id)"
            ")"
        );
        Execute(
            "create table if not exists environment_fingerprint ("
            "id ubigint primary key,"
            "smarts varchar not null,"
            "pseudosmiles varchar not null,"
            "parent_smarts varchar not null,"
            "unique (smarts, pseudosmiles, parent_smarts)"
            ")"
        );
        Execute(
            "create table if not exists rule_environment ("
            "id ubigint primary key,"
            "rule_id ubigint not null references rule(id),"
            "environment_fingerprint_id ubigint not null references environment_fingerprint(id),"
            "radius integer not null,"
            "num_pairs uinteger not null,"
            "unique (rule_id, environment_fingerprint_id, radius)"
            ")"
        );
        Execute(
            "create table if not exists constant_smiles ("
            "id ubigint primary key,"
            "smiles varchar not null unique"
            ")"
        );
        Execute(
            "create table if not exists pair ("
            "id ubigint primary key,"
            "rule_environment_id ubigint not null references rule_environment(id),"
            "compound1_id ubigint not null references compound(id),"
            "compound2_id ubigint not null references compound(id),"
            "constant_id ubigint not null references constant_smiles(id),"
            "cut_count uinteger not null,"
            "heavy_atom_delta integer not null,"
            "heavy_bond_delta integer not null"
            ")"
        );
        Execute("commit");
    } catch (...) {
        try {
            Execute("rollback");
        } catch (const StorageError&) {
        }
        throw;
    }
}

void DuckDBStore::Execute(const std::string& sql) {
    if (!connection_) {
        throw StorageError("DuckDB connection is not open");
    }

    std::unique_ptr<duckdb::QueryResult> result = connection_->Query(sql);
    if (!result) {
        throw StorageError("DuckDB query returned no result: " + sql);
    }
    require_success(*result, sql);
}

void DuckDBStore::AddMolecule(const MoleculeRecord& molecule) {
    if (molecule.GetCanonicalSmiles().empty()) {
        throw StorageError("molecule record has no canonical SMILES");
    }
    if (HasMolecule(molecule.GetInternalId())) {
        throw StorageError(
            "duplicate internal molecule id: " + std::to_string(molecule.GetInternalId())
        );
    }
    if (
        molecule.HasExternalId() &&
        has_external_molecule_id(connection_, molecule.GetExternalId())
    ) {
        throw DuplicateIdError("duplicate molecule id: " + molecule.GetExternalId());
    }

    const std::string sql =
        "insert into compound ("
        "id, public_id, input_smiles, clean_smiles, clean_num_heavies"
        ") values (?, ?, ?, ?, ?)";
    duckdb::vector<duckdb::Value> values = {
        duckdb::Value::UBIGINT(molecule.GetInternalId()),
        string_or_null(molecule.GetExternalId()),
        duckdb::Value(molecule.GetCanonicalSmiles()),
        duckdb::Value(molecule.GetCanonicalSmiles()),
        duckdb::Value::UINTEGER(molecule.GetHeavyAtomCount()),
    };
    execute_prepared(connection_, sql, std::move(values));
}

LoadReport DuckDBStore::AddMoleculesFromSmilesFile(const std::string& smiles_path) {
    InitializeSchema();

    std::ifstream input(smiles_path);
    if (!input) {
        throw StorageError("failed to open SMILES file: " + smiles_path);
    }

    LoadReport report;
    std::uint64_t next_id = get_next_id(connection_, "compound", "id");
    std::string line;
    unsigned int row_number = 0;
    while (std::getline(input, line)) {
        ++row_number;
        const std::string stripped = trim_copy(line);
        if (stripped.empty() || stripped[0] == '#') {
            continue;
        }

        std::istringstream row_stream(stripped);
        std::string smiles;
        std::string external_id;
        row_stream >> smiles;
        if (!(row_stream >> external_id)) {
            external_id = make_generated_external_id(next_id);
        }

        try {
            const MoleculeRecord molecule = MoleculeRecord::FromSmiles(
                static_cast<unsigned int>(next_id),
                smiles,
                external_id
            );
            AddMolecule(molecule);
        } catch (const std::exception& exc) {
            report.RecordRejected(row_number, exc.what());
            continue;
        }

        report.RecordAccepted(external_id);
        ++next_id;
    }

    return report;
}

LoadReport DuckDBStore::AddPropertiesFromCsvFile(
    const std::string& csv_path,
    const std::string& id_column,
    const std::vector<std::string>& property_columns
) {
    InitializeSchema();

    std::ifstream input(csv_path);
    if (!input) {
        throw StorageError("failed to open property CSV file: " + csv_path);
    }

    std::string line;
    if (!std::getline(input, line)) {
        throw StorageError("property CSV file must contain a header row: " + csv_path);
    }

    const std::vector<std::string> header = parse_csv_line(line);
    const std::unordered_map<std::string, std::size_t> header_index =
        build_header_index(header);
    const std::string resolved_id_column = resolve_id_column(header_index, id_column);
    const std::vector<std::string> resolved_property_columns =
        resolve_property_columns(
            header,
            header_index,
            resolved_id_column,
            property_columns
        );
    const std::size_t id_index = header_index.at(resolved_id_column);

    LoadReport report;
    unsigned int row_number = 1;
    while (std::getline(input, line)) {
        ++row_number;
        if (trim_copy(line).empty()) {
            continue;
        }

        try {
            const std::vector<std::string> fields = parse_csv_line(line);
            if (fields.size() != header.size()) {
                throw StorageError(
                    "CSV row has " + std::to_string(fields.size()) +
                    " fields but header has " + std::to_string(header.size())
                );
            }

            const std::string external_id = fields[id_index];
            if (external_id.empty()) {
                throw StorageError("CSV id value must not be blank");
            }

            const std::uint64_t molecule_id =
                find_molecule_internal_id_by_external_id(connection_, external_id);
            if (molecule_id == 0) {
                throw StorageError("unknown molecule id in property CSV: " + external_id);
            }

            std::vector<std::pair<std::string, double>> parsed_values;
            for (const std::string& property_name : resolved_property_columns) {
                const std::string value_text = fields[header_index.at(property_name)];
                if (value_text.empty() || value_text == "*") {
                    continue;
                }
                parsed_values.push_back({
                    property_name,
                    parse_property_value(property_name, value_text),
                });
            }

            for (const auto& property_value : parsed_values) {
                AddMoleculeProperty(
                    static_cast<unsigned int>(molecule_id),
                    property_value.first,
                    property_value.second
                );
            }
            report.RecordAccepted(external_id);
        } catch (const std::exception& exc) {
            report.RecordRejected(row_number, exc.what());
        }
    }

    return report;
}

LoadReport DuckDBStore::AddPropertiesFromCsvFile(
    const std::string& csv_path,
    const std::string& id_column
) {
    return AddPropertiesFromCsvFile(csv_path, id_column, std::vector<std::string>());
}

LoadReport DuckDBStore::AddPropertiesFromCsvFile(const std::string& csv_path) {
    return AddPropertiesFromCsvFile(csv_path, "id", std::vector<std::string>());
}

void DuckDBStore::AddMoleculeProperty(
    unsigned int molecule_id,
    const std::string& property_name,
    double value
) {
    if (property_name.empty()) {
        throw StorageError("molecule property name must not be empty");
    }
    if (!HasMolecule(molecule_id)) {
        throw StorageError("cannot add property for unknown molecule");
    }

    const std::uint64_t property_name_id = get_or_create_named_row_id(
        connection_,
        "property_name",
        "id",
        "name",
        property_name
    );
    const std::uint64_t property_id = get_next_id(connection_, "compound_property", "id");
    const std::string sql =
        "insert into compound_property (id, compound_id, property_name_id, value) "
        "values (?, ?, ?, ?) "
        "on conflict (compound_id, property_name_id) do update "
        "set value = excluded.value";
    duckdb::vector<duckdb::Value> values = {
        duckdb::Value::UBIGINT(property_id),
        duckdb::Value::UBIGINT(molecule_id),
        duckdb::Value::UBIGINT(property_name_id),
        duckdb::Value::DOUBLE(value),
    };
    execute_prepared(connection_, sql, std::move(values));
}

void DuckDBStore::AddPair(const MatchedPair& pair) {
    if (!HasMolecule(pair.GetSourceMoleculeId())) {
        throw StorageError("cannot add pair with unknown source molecule");
    }
    if (!HasMolecule(pair.GetTargetMoleculeId())) {
        throw StorageError("cannot add pair with unknown target molecule");
    }

    const std::uint64_t constant_id = get_or_create_named_row_id(
        connection_,
        "constant_smiles",
        "id",
        "smiles",
        pair.GetConstantSmiles()
    );
    const std::uint64_t from_smiles_id = get_or_create_named_row_id(
        connection_,
        "rule_smiles",
        "id",
        "smiles",
        pair.GetSourceVariableSmiles()
    );
    const std::uint64_t to_smiles_id = get_or_create_named_row_id(
        connection_,
        "rule_smiles",
        "id",
        "smiles",
        pair.GetTargetVariableSmiles()
    );
    const std::uint64_t rule_id = get_or_create_rule_id(
        connection_,
        from_smiles_id,
        to_smiles_id
    );
    // Radius-zero empty environments preserve the MMPDB rule-environment
    // boundary until atom-context fingerprinting is implemented.
    const std::uint64_t environment_fingerprint_id =
        get_or_create_environment_fingerprint_id(connection_, "", "", "");
    const std::uint64_t rule_environment_id =
        get_or_create_rule_environment_id(connection_, rule_id, environment_fingerprint_id, 0);
    const std::uint64_t pair_id = get_next_id(connection_, "pair", "id");

    const std::string sql =
        "insert into pair ("
        "id, rule_environment_id, compound1_id, compound2_id, constant_id, "
        "cut_count, heavy_atom_delta, heavy_bond_delta"
        ") values (?, ?, ?, ?, ?, ?, ?, ?)";
    duckdb::vector<duckdb::Value> values = {
        duckdb::Value::UBIGINT(pair_id),
        duckdb::Value::UBIGINT(rule_environment_id),
        duckdb::Value::UBIGINT(pair.GetSourceMoleculeId()),
        duckdb::Value::UBIGINT(pair.GetTargetMoleculeId()),
        duckdb::Value::UBIGINT(constant_id),
        duckdb::Value::UINTEGER(pair.GetCutCount()),
        duckdb::Value::INTEGER(pair.GetHeavyAtomDelta()),
        duckdb::Value::INTEGER(pair.GetHeavyBondDelta()),
    };
    execute_prepared(connection_, sql, std::move(values));
    increment_rule_environment_pair_count(connection_, rule_environment_id);
}

void DuckDBStore::AddPairs(const std::vector<MatchedPair>& pairs) {
    Execute("begin transaction");
    try {
        for (const MatchedPair& pair : pairs) {
            AddPair(pair);
        }
        Execute("commit");
    } catch (...) {
        try {
            Execute("rollback");
        } catch (const StorageError&) {
        }
        throw;
    }
}

bool DuckDBStore::HasTable(const std::string& table_name) const {
    const std::vector<std::string> tables = GetTableNames();
    return std::find(tables.begin(), tables.end(), table_name) != tables.end();
}

bool DuckDBStore::HasMolecule(unsigned int internal_id) const {
    const std::string sql = "select count(*) from compound where id = ?";
    duckdb::vector<duckdb::Value> values = {duckdb::Value::UBIGINT(internal_id)};
    std::unique_ptr<duckdb::QueryResult> result =
        execute_prepared(connection_, sql, std::move(values));

    for (const auto& row : *result) {
        return row.GetValue<std::int64_t>(0) > 0;
    }
    return false;
}

std::uint64_t DuckDBStore::GetRowCount(const std::string& table_name) const {
    if (!is_base_table_name(table_name)) {
        throw StorageError("unknown DuckDB base table: " + table_name);
    }

    const std::string sql = "select count(*) from " + resolve_table_alias(table_name);
    std::unique_ptr<duckdb::QueryResult> result = connection_->Query(sql);
    if (!result) {
        throw StorageError("DuckDB query returned no result: " + sql);
    }
    require_success(*result, sql);

    for (const auto& row : *result) {
        return static_cast<std::uint64_t>(row.GetValue<std::int64_t>(0));
    }
    return 0;
}

double DuckDBStore::GetMoleculeProperty(
    unsigned int molecule_id,
    const std::string& property_name
) const {
    if (property_name.empty()) {
        throw StorageError("molecule property name must not be empty");
    }

    const std::string sql =
        "select compound_property.value "
        "from compound_property "
        "join property_name on property_name.id = compound_property.property_name_id "
        "where compound_property.compound_id = ? and property_name.name = ?";
    duckdb::vector<duckdb::Value> values = {
        duckdb::Value::UBIGINT(molecule_id),
        duckdb::Value(property_name),
    };
    std::unique_ptr<duckdb::QueryResult> result =
        execute_prepared(connection_, sql, std::move(values));

    for (const auto& row : *result) {
        return row.GetValue<double>(0);
    }

    throw StorageError("molecule property not found");
}

std::vector<MatchedPair> DuckDBStore::GetPairs() const {
    return GetPairs(QueryOptions());
}

std::vector<MatchedPair> DuckDBStore::GetPairs(const QueryOptions& options) const {
    const std::string sql = build_pair_query(options);
    std::unique_ptr<duckdb::QueryResult> result = connection_->Query(sql);
    if (!result) {
        throw StorageError("DuckDB query returned no result: " + sql);
    }
    require_success(*result, sql);

    std::vector<MatchedPair> pairs;
    for (const auto& row : *result) {
        MatchedPair pair(
            static_cast<unsigned int>(row.GetValue<std::uint64_t>(0)),
            static_cast<unsigned int>(row.GetValue<std::uint64_t>(1)),
            row.GetValue<std::string>(2),
            row.GetValue<std::string>(3),
            row.GetValue<std::string>(4),
            row.GetValue<std::string>(5),
            row.GetValue<std::string>(6),
            row.GetValue<std::string>(7),
            row.GetValue<std::string>(8),
            static_cast<unsigned int>(row.GetValue<std::uint32_t>(9)),
            row.GetValue<std::int32_t>(10),
            row.GetValue<std::int32_t>(11)
        );

        const std::string property_sql =
            "select property_name.name, source_property.value, target_property.value "
            "from compound_property source_property "
            "join compound_property target_property "
            "on target_property.property_name_id = source_property.property_name_id "
            "join property_name on property_name.id = source_property.property_name_id "
            "where source_property.compound_id = ? and target_property.compound_id = ? "
            "order by property_name.name";
        duckdb::vector<duckdb::Value> property_values = {
            duckdb::Value::UBIGINT(pair.GetSourceMoleculeId()),
            duckdb::Value::UBIGINT(pair.GetTargetMoleculeId()),
        };
        std::unique_ptr<duckdb::QueryResult> property_result =
            execute_prepared(connection_, property_sql, std::move(property_values));
        for (const auto& property_row : *property_result) {
            pair.SetProperty(
                property_row.GetValue<std::string>(0),
                property_row.GetValue<double>(1),
                property_row.GetValue<double>(2)
            );
        }

        pairs.push_back(pair);
    }
    return apply_pair_scoring(pairs, options.GetScoringOptions());
}

std::vector<Transform> DuckDBStore::GetTransforms() const {
    return GetTransforms(QueryOptions());
}

std::vector<Transform> DuckDBStore::GetTransforms(const QueryOptions& options) const {
    std::map<std::string, Transform> transforms_by_smiles;
    for (const MatchedPair& pair : GetPairs(options)) {
        const std::string& transform_smiles = pair.GetTransformSmiles();
        auto inserted = transforms_by_smiles.emplace(
            transform_smiles,
            Transform(transform_smiles)
        );
        inserted.first->second.AddPair(pair);
    }

    std::vector<Transform> transforms;
    transforms.reserve(transforms_by_smiles.size());
    for (const auto& entry : transforms_by_smiles) {
        transforms.push_back(entry.second);
    }
    return transforms;
}

std::vector<std::string> DuckDBStore::GetTableNames() const {
    if (!connection_) {
        throw StorageError("DuckDB connection is not open");
    }

    const std::string sql =
        "select table_name "
        "from information_schema.tables "
        "where table_schema = 'main' and table_type = 'BASE TABLE' "
        "order by table_name";
    std::unique_ptr<duckdb::QueryResult> result = connection_->Query(sql);
    if (!result) {
        throw StorageError("DuckDB query returned no result: " + sql);
    }
    require_success(*result, sql);

    std::vector<std::string> table_names;
    for (const auto& row : *result) {
        table_names.push_back(row.GetValue<std::string>(0));
    }
    return table_names;
}

}  // namespace OEMMPA
