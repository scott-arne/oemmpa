#include "oemmpa/DuckDBStore.h"

#include "oedesalt/Desalter.h"
#include "oemmpa/EnvironmentFingerprint.h"
#include "oemmpa/Error.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/PairScoring.h"

#include <duckdb.hpp>

#include <oechem.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
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
        "constant_environment",
        "constant_smiles",
        "dataset",
        "environment_fingerprint",
        "pair",
        "property_name",
        "rule",
        "rule_environment",
        "rule_environment_statistics",
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

// Staged new-row payloads for the bulk resolve-then-append path. They live in
// the anonymous namespace so the file-local Append* helpers below can name
// them; only rows not already present in the store are staged for append.
struct NewConstant {
    std::uint64_t id;
    std::string smiles;
};

struct NewRuleSmiles {
    std::uint64_t id;
    std::string smiles;
};

struct NewRule {
    std::uint64_t id;
    std::uint64_t from_id;
    std::uint64_t to_id;
};

struct NewFingerprint {
    std::uint64_t id;
    std::string smarts;
    std::string pseudo;
    std::string parent;
};

struct NewRuleEnvironment {
    std::uint64_t id;
    std::uint64_t rule_id;
    std::uint64_t fingerprint_id;
    int radius;
    std::uint64_t num_pairs;
};

struct PairRow {
    std::uint64_t id;
    std::uint64_t rule_id;
    std::uint64_t constant_id;
    std::uint64_t compound1_id;
    std::uint64_t compound2_id;
    unsigned int cut_count;
    int heavy_atom_delta;
    int heavy_bond_delta;
};

struct NewConstantEnvironment {
    std::uint64_t constant_id;
    int radius;
    std::uint64_t fingerprint_id;
};

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

// Name of the id-allocation sequence backing a table's primary key.
std::string id_sequence_name(const std::string& table_name) {
    return "seq_" + table_name + "_id";
}

// Current maximum id in a table, or 0 when empty. Used once per table at
// schema-init time to seed its id sequence; not on the per-insert hot path.
std::uint64_t get_max_id(
    const std::unique_ptr<duckdb::Connection>& connection,
    const std::string& table_name
) {
    const std::string sql = "select coalesce(max(id), 0) from " + table_name;
    std::unique_ptr<duckdb::QueryResult> result = connection->Query(sql);
    if (!result) {
        throw StorageError("DuckDB query returned no result: " + sql);
    }
    require_success(*result, sql);

    for (const auto& row : *result) {
        return row.GetValue<std::uint64_t>(0);
    }
    throw StorageError("DuckDB did not return a max ID for " + table_name);
}

// Allocate the next primary-key id for a table. Backed by a DuckDB sequence
// (see InitializeSchema), so this is an atomic O(1) counter rather than a
// repeated max(id) table scan.
std::uint64_t get_next_id(
    const std::unique_ptr<duckdb::Connection>& connection,
    const std::string& table_name,
    const std::string& /*id_column*/
) {
    const std::string sql = "select nextval('" + id_sequence_name(table_name) + "')";
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


struct AggregateStatistics {
    std::uint32_t count = 0;
    double avg = 0.0;
    duckdb::Value std = duckdb::Value(nullptr);
    duckdb::Value kurtosis = duckdb::Value(nullptr);
    duckdb::Value skewness = duckdb::Value(nullptr);
    double min = 0.0;
    double q1 = 0.0;
    double median = 0.0;
    double q3 = 0.0;
    double max = 0.0;
    duckdb::Value paired_t = duckdb::Value(nullptr);
    duckdb::Value p_value = duckdb::Value(nullptr);
};

// Continued-fraction evaluation for the regularized incomplete beta function
// (Numerical Recipes betacf, modified Lentz). Used to compute the Student's t
// two-sided p-value without a third-party dependency, matching scipy's
// stats.t.sf(|t|, df) * 2 to machine precision.
double incomplete_beta_cf(double a, double b, double x) {
    const int MAX_ITERATIONS = 200;
    const double EPSILON = 3.0e-16;
    const double FLOOR = 1.0e-300;
    const double qab = a + b;
    const double qap = a + 1.0;
    const double qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::fabs(d) < FLOOR) {
        d = FLOOR;
    }
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m <= MAX_ITERATIONS; ++m) {
        const double m2 = 2.0 * m;
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::fabs(d) < FLOOR) {
            d = FLOOR;
        }
        c = 1.0 + aa / c;
        if (std::fabs(c) < FLOOR) {
            c = FLOOR;
        }
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (std::fabs(d) < FLOOR) {
            d = FLOOR;
        }
        c = 1.0 + aa / c;
        if (std::fabs(c) < FLOOR) {
            c = FLOOR;
        }
        d = 1.0 / d;
        const double delta = d * c;
        h *= delta;
        if (std::fabs(delta - 1.0) < EPSILON) {
            break;
        }
    }
    return h;
}

// Regularized incomplete beta function I_x(a, b).
double regularized_incomplete_beta(double a, double b, double x) {
    if (x <= 0.0) {
        return 0.0;
    }
    if (x >= 1.0) {
        return 1.0;
    }
    const double front = std::exp(
        std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) +
        a * std::log(x) + b * std::log1p(-x)
    );
    if (x < (a + 1.0) / (a + b + 2.0)) {
        return front * incomplete_beta_cf(a, b, x) / a;
    }
    return 1.0 - front * incomplete_beta_cf(b, a, 1.0 - x) / b;
}

// Two-sided Student's t p-value for a t statistic with the given degrees of
// freedom. Equals scipy.stats.t.sf(|t|, df) * 2 via the identity
// 2 * sf(|t|, df) == I_{df/(df+t^2)}(df/2, 1/2).
double two_sided_t_p_value(double t, double df) {
    if (df <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double x = df / (df + t * t);
    return regularized_incomplete_beta(df / 2.0, 0.5, x);
}

double median(const std::vector<double>& values) {
    const std::size_t count = values.size();
    if (count == 0) {
        throw StorageError("cannot compute median for an empty value set");
    }

    const std::size_t half = count / 2;
    if (count % 2 == 1) {
        return values[half];
    }
    return (values[half - 1] + values[half]) / 2.0;
}

double median_range(
    const std::vector<double>& values,
    std::size_t begin,
    std::size_t end
) {
    std::vector<double> subset(values.begin() + begin, values.begin() + end);
    return median(subset);
}

std::tuple<double, double, double> quartiles(const std::vector<double>& values) {
    const std::size_t count = values.size();
    if (count == 1) {
        return {values[0], values[0], values[0]};
    }

    const double median_value = median(values);
    const std::size_t half = count / 2;
    if (count % 2 == 0) {
        return {
            median_range(values, 0, half),
            median_value,
            median_range(values, half, count),
        };
    }

    if (count % 4 == 1) {
        const std::size_t middle = (count - 1) / 4;
        const double q1 = 0.25 * values[middle - 1] + 0.75 * values[middle];
        const double q3 =
            0.75 * values[3 * middle] + 0.25 * values[3 * middle + 1];
        return {q1, median_value, q3};
    }

    const std::size_t middle = (count - 3) / 4;
    const double q1 = 0.75 * values[middle] + 0.25 * values[middle + 1];
    const double q3 =
        0.25 * values[3 * middle + 1] + 0.75 * values[3 * middle + 2];
    return {q1, median_value, q3};
}

duckdb::Value sample_standard_deviation(const std::vector<double>& values) {
    std::uint32_t count = 0;
    double mean = 0.0;
    double m2 = 0.0;
    for (const double value : values) {
        ++count;
        const double delta = value - mean;
        mean += delta / static_cast<double>(count);
        m2 += delta * (value - mean);
    }
    if (count < 2) {
        return duckdb::Value(nullptr);
    }
    return duckdb::Value::DOUBLE(std::sqrt(m2 / static_cast<double>(count - 1)));
}

duckdb::Value kurtosis(const std::vector<double>& values) {
    std::uint32_t count = 0;
    double mean = 0.0;
    double m2 = 0.0;
    double m3 = 0.0;
    double m4 = 0.0;
    for (const double value : values) {
        const std::uint32_t previous_count = count;
        ++count;
        const double delta = value - mean;
        const double delta_n = delta / static_cast<double>(count);
        const double delta_n2 = delta_n * delta_n;
        const double term1 = delta * delta_n * static_cast<double>(previous_count);
        mean += delta_n;
        m4 +=
            term1 * delta_n2 *
                static_cast<double>(count * count - 3 * count + 3) +
            6.0 * delta_n2 * m2 -
            4.0 * delta_n * m3;
        m3 +=
            term1 * delta_n * static_cast<double>(count - 2) -
            3.0 * delta_n * m2;
        m2 += term1;
    }
    if (m2 == 0.0) {
        return duckdb::Value(nullptr);
    }
    return duckdb::Value::DOUBLE(
        (static_cast<double>(count) * m4) / (m2 * m2) - 3.0
    );
}

duckdb::Value skewness(const std::vector<double>& values, double avg) {
    const std::size_t count = values.size();
    double skew_top = 0.0;
    double squared_delta_sum = 0.0;
    for (const double value : values) {
        const double delta = value - avg;
        skew_top += delta * delta * delta;
        squared_delta_sum += delta * delta;
    }

    skew_top /= static_cast<double>(count);
    if (skew_top == 0.0) {
        return duckdb::Value::DOUBLE(0.0);
    }

    const double skew_bot = std::pow(
        squared_delta_sum / static_cast<double>(count - 1),
        1.5
    );
    if (skew_bot == 0.0) {
        return duckdb::Value(nullptr);
    }
    return duckdb::Value::DOUBLE(skew_top / skew_bot);
}

AggregateStatistics aggregate_values(std::vector<double> values) {
    if (values.empty()) {
        throw StorageError("cannot compute statistics for an empty value set");
    }

    std::sort(values.begin(), values.end());

    AggregateStatistics statistics;
    statistics.count = static_cast<std::uint32_t>(values.size());
    double sum = 0.0;
    for (const double value : values) {
        sum += value;
    }
    statistics.avg = sum / static_cast<double>(statistics.count);
    statistics.std = sample_standard_deviation(values);
    if (statistics.count > 2) {
        statistics.kurtosis = kurtosis(values);
        statistics.skewness = skewness(values, statistics.avg);
    }
    std::tie(statistics.q1, statistics.median, statistics.q3) = quartiles(values);
    statistics.min = values.front();
    statistics.max = values.back();

    if (statistics.count > 1) {
        const double std_value = statistics.std.GetValue<double>();
        if (std_value == 0.0) {
            // Zero variance: paired_t saturates and the p-value is undefined,
            // so it is left NULL (matching the Python _p_value contract).
            statistics.paired_t = duckdb::Value::DOUBLE(100000000.0);
        } else {
            const double paired_t = std::min(
                (statistics.avg / std_value) *
                    std::sqrt(static_cast<double>(statistics.count)),
                100000000.0
            );
            statistics.paired_t = duckdb::Value::DOUBLE(paired_t);
            // Two-sided p-value from the (clamped) t statistic with
            // df = count - 1, matching Python's scipy-based computation.
            statistics.p_value = duckdb::Value::DOUBLE(
                two_sided_t_p_value(
                    paired_t,
                    static_cast<double>(statistics.count - 1)
                )
            );
        }
    }

    return statistics;
}

using RuleEnvironmentPropertyKey = std::pair<std::uint64_t, std::uint64_t>;

std::map<RuleEnvironmentPropertyKey, std::vector<double>>
collect_rule_environment_property_deltas(
    const std::unique_ptr<duckdb::Connection>& connection
) {
    const std::string sql =
        "select "
        "rule_env.id as rule_environment_id, "
        "source_property.property_name_id, "
        "target_property.value - source_property.value "
        "from pair "
        "join constant_environment ce on ce.constant_id = pair.constant_id "
        "join rule_environment rule_env "
        "on rule_env.rule_id = pair.rule_id "
        "and rule_env.environment_fingerprint_id = ce.environment_fingerprint_id "
        "and rule_env.radius = ce.radius "
        "join compound_property source_property "
        "on source_property.compound_id = pair.compound1_id "
        "join compound_property target_property "
        "on target_property.compound_id = pair.compound2_id "
        "and target_property.property_name_id = source_property.property_name_id "
        "order by rule_env.id, source_property.property_name_id, pair.id";
    std::unique_ptr<duckdb::QueryResult> result = connection->Query(sql);
    if (!result) {
        throw StorageError("DuckDB query returned no result: " + sql);
    }
    require_success(*result, sql);

    std::map<RuleEnvironmentPropertyKey, std::vector<double>> deltas_by_key;
    for (const auto& row : *result) {
        deltas_by_key[{
            row.GetValue<std::uint64_t>(0),
            row.GetValue<std::uint64_t>(1),
        }].push_back(row.GetValue<double>(2));
    }
    return deltas_by_key;
}

void refresh_dataset_counts(
    const std::unique_ptr<duckdb::Connection>& connection
) {
    const std::string insert_sql =
        "insert into dataset ("
        "id, oemmpa_schema_version, title, fragment_options, index_options, is_symmetric"
        ") values (1, " + std::to_string(DuckDBStore::kOemmpaSchemaVersion) +
        ", '', '', '', true) "
        "on conflict (id) do nothing";
    std::unique_ptr<duckdb::QueryResult> insert_result = connection->Query(insert_sql);
    if (!insert_result) {
        throw StorageError("DuckDB query returned no result: " + insert_sql);
    }
    require_success(*insert_result, insert_sql);

    const std::string update_sql =
        "update dataset set "
        "num_compounds = (select count(*) from compound),"
        "num_rules = (select count(*) from rule),"
        "num_pairs = (select count(*) from pair),"
        "num_rule_environments = (select count(*) from rule_environment),"
        "num_rule_environment_stats = (select count(*) from rule_environment_statistics) "
        "where id = 1";
    std::unique_ptr<duckdb::QueryResult> update_result = connection->Query(update_sql);
    if (!update_result) {
        throw StorageError("DuckDB query returned no result: " + update_sql);
    }
    require_success(*update_result, update_sql);
}

void refresh_rule_environment_statistics(
    const std::unique_ptr<duckdb::Connection>& connection
) {
    const std::map<RuleEnvironmentPropertyKey, std::vector<double>> deltas_by_key =
        collect_rule_environment_property_deltas(connection);

    const std::string delete_sql = "delete from rule_environment_statistics";
    std::unique_ptr<duckdb::QueryResult> delete_result = connection->Query(delete_sql);
    if (!delete_result) {
        throw StorageError("DuckDB query returned no result: " + delete_sql);
    }
    require_success(*delete_result, delete_sql);

    const std::string insert_sql =
        "insert into rule_environment_statistics ("
        "id, rule_environment_id, property_name_id, count, avg, std, kurtosis, skewness, "
        "min, q1, median, q3, max, paired_t, p_value"
        ") values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";
    std::uint64_t next_id = 1;
    for (const auto& entry : deltas_by_key) {
        const AggregateStatistics statistics = aggregate_values(entry.second);
        duckdb::vector<duckdb::Value> values = {
            duckdb::Value::UBIGINT(next_id),
            duckdb::Value::UBIGINT(entry.first.first),
            duckdb::Value::UBIGINT(entry.first.second),
            duckdb::Value::UINTEGER(statistics.count),
            duckdb::Value::DOUBLE(statistics.avg),
            statistics.std,
            statistics.kurtosis,
            statistics.skewness,
            duckdb::Value::DOUBLE(statistics.min),
            duckdb::Value::DOUBLE(statistics.q1),
            duckdb::Value::DOUBLE(statistics.median),
            duckdb::Value::DOUBLE(statistics.q3),
            duckdb::Value::DOUBLE(statistics.max),
            statistics.paired_t,
            statistics.p_value,
        };
        execute_prepared(connection, insert_sql, std::move(values));
        ++next_id;
    }
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

// The synthesized hydrogen pseudo-fragment. MMPDB (and the in-memory query
// path) exempt this fragment from variable-size bounds, so the store must too.
const char* const kHydrogenVariableSmiles = "[*:1][H]";

// Heavy-atom count of a variable SMILES, counted the same way as the in-memory
// path (OEChem OEIsHeavy). Cached by SMILES so the common case of many pairs
// sharing a variable fragment parses each distinct fragment once.
unsigned int variable_heavy_atom_count(
    const std::string& variable_smiles,
    std::unordered_map<std::string, unsigned int>& cache
) {
    const auto cached = cache.find(variable_smiles);
    if (cached != cache.end()) {
        return cached->second;
    }
    OEChem::OEGraphMol mol;
    if (!OEChem::OESmilesToMol(mol, variable_smiles)) {
        throw StorageError("invalid variable SMILES in store: " + variable_smiles);
    }
    const unsigned int count = OEChem::OECount(mol, OEChem::OEIsHeavy());
    cache.emplace(variable_smiles, count);
    return count;
}

// Apply the variable-fragment size bounds to already-materialized store pairs,
// using the SAME QueryOptions::AllowsVariableFragment predicate as the
// in-memory MemoryIndex path so both backends filter identically. |V| is
// computed from the pair's variable SMILES; the ratio denominator is the whole
// molecule's heavy count, fetched for the referenced compounds in one query.
// The [*:1][H] pseudo-fragment is exempt per side (matching MMPDB), so a min
// bound never drops an H<->heavy substitution.
std::vector<MatchedPair> filter_pairs_by_variable_bounds(
    const std::unique_ptr<duckdb::Connection>& connection,
    std::vector<MatchedPair> pairs,
    const QueryOptions& options
) {
    if (!options.HasVariableFragmentBounds() || pairs.empty()) {
        return pairs;
    }

    // Fetch whole-molecule heavy counts for every referenced compound in one
    // query (mirrors read_pairs_with_properties' batched-lookup pattern).
    std::set<std::uint64_t> compound_ids;
    for (const MatchedPair& pair : pairs) {
        compound_ids.insert(pair.GetSourceMoleculeId());
        compound_ids.insert(pair.GetTargetMoleculeId());
    }
    std::ostringstream id_list;
    bool first = true;
    for (const std::uint64_t compound_id : compound_ids) {
        if (!first) {
            id_list << ",";
        }
        id_list << compound_id;
        first = false;
    }
    const std::string heavies_sql =
        "select id, clean_num_heavies from compound where id in (" +
        id_list.str() + ")";
    std::unique_ptr<duckdb::QueryResult> heavies_result =
        connection->Query(heavies_sql);
    if (!heavies_result) {
        throw StorageError("DuckDB query returned no result: " + heavies_sql);
    }
    require_success(*heavies_result, heavies_sql);
    std::unordered_map<std::uint64_t, unsigned int> heavies_by_compound;
    for (const auto& row : *heavies_result) {
        heavies_by_compound.emplace(
            row.GetValue<std::uint64_t>(0),
            static_cast<unsigned int>(row.GetValue<std::uint32_t>(1)));
    }

    std::unordered_map<std::string, unsigned int> variable_heavy_cache;
    const auto side_allowed = [&](const std::string& variable_smiles,
                                  std::uint64_t molecule_id) {
        if (variable_smiles == kHydrogenVariableSmiles) {
            return true;  // pseudo-fragment is exempt, matching MMPDB
        }
        const auto heavies = heavies_by_compound.find(molecule_id);
        const unsigned int molecule_heavies =
            heavies == heavies_by_compound.end() ? 0U : heavies->second;
        return options.AllowsVariableFragment(
            variable_heavy_atom_count(variable_smiles, variable_heavy_cache),
            molecule_heavies);
    };

    std::vector<MatchedPair> kept;
    kept.reserve(pairs.size());
    for (MatchedPair& pair : pairs) {
        if (side_allowed(pair.GetSourceVariableSmiles(), pair.GetSourceMoleculeId()) &&
            side_allowed(pair.GetTargetVariableSmiles(), pair.GetTargetMoleculeId())) {
            kept.push_back(std::move(pair));
        }
    }
    return kept;
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

std::string build_pair_query(
    const QueryOptions& options,
    std::uint64_t rule_environment_id = 0
) {
    std::ostringstream sql;
    sql << "select "
        << "p.compound1_id, p.compound2_id, "
        << "coalesce(source_molecule.public_id, '') as source_public_id, "
        << "coalesce(target_molecule.public_id, '') as target_public_id, "
        << "source_molecule.clean_smiles as source_smiles, "
        << "target_molecule.clean_smiles as target_smiles, "
        << "c.smiles as constant_smiles, "
        << "source_variable.smiles as source_variable_smiles, "
        << "target_variable.smiles as target_variable_smiles, "
        << "p.cut_count, p.heavy_atom_delta, p.heavy_bond_delta "
        << "from pair p "
        << "join compound source_molecule on source_molecule.id = p.compound1_id "
        << "join compound target_molecule on target_molecule.id = p.compound2_id "
        << "join constant_smiles c on c.id = p.constant_id "
        << "join rule r on r.id = p.rule_id "
        << "join rule_smiles source_variable on source_variable.id = r.from_smiles_id "
        << "join rule_smiles target_variable on target_variable.id = r.to_smiles_id ";

    if (rule_environment_id > 0) {
        // Restrict to physical pairs whose (rule, constant) resolve to this
        // environment's fingerprint at its radius.
        sql << "join constant_environment ce on ce.constant_id = p.constant_id "
            << "join rule_environment rule_env "
            << "on rule_env.rule_id = p.rule_id "
            << "and rule_env.environment_fingerprint_id = ce.environment_fingerprint_id "
            << "and rule_env.radius = ce.radius ";
    }

    sql << "where true ";
    if (rule_environment_id > 0) {
        sql << "and rule_env.id = " << rule_environment_id << " ";
    }
    if (!options.GetSymmetric()) {
        sql << "and p.compound1_id < p.compound2_id ";
    }
    if (options.GetMaxHeavyAtomChange() >= 0) {
        sql << "and abs(p.heavy_atom_delta) <= "
            << options.GetMaxHeavyAtomChange() << " ";
    }
    if (options.GetMaxRelativeHeavyAtomChange() >= 0.0) {
        sql << "and case "
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

std::string build_rule_environment_statistics_query(bool filter_property) {
    std::ostringstream sql;
    sql
        << "select "
        << "stats.rule_environment_id, "
        << "property_name.name, "
        << "from_smiles.smiles, "
        << "to_smiles.smiles, "
        << "rule_env.radius, "
        << "coalesce(environment_fingerprint.smarts, ''), "
        << "coalesce(environment_fingerprint.pseudosmiles, ''), "
        << "coalesce(environment_fingerprint.parent_smarts, ''), "
        << "stats.count, "
        << "stats.avg, "
        << "stats.std is not null, coalesce(stats.std, 0.0), "
        << "stats.kurtosis is not null, coalesce(stats.kurtosis, 0.0), "
        << "stats.skewness is not null, coalesce(stats.skewness, 0.0), "
        << "stats.min, stats.q1, stats.median, stats.q3, stats.max, "
        << "stats.paired_t is not null, coalesce(stats.paired_t, 0.0), "
        << "stats.p_value is not null, coalesce(stats.p_value, 0.0) "
        << "from rule_environment_statistics stats "
        << "join property_name on property_name.id = stats.property_name_id "
        << "join rule_environment rule_env on rule_env.id = stats.rule_environment_id "
        << "join environment_fingerprint "
        << "on environment_fingerprint.id = rule_env.environment_fingerprint_id "
        << "join rule on rule.id = rule_env.rule_id "
        << "join rule_smiles from_smiles on from_smiles.id = rule.from_smiles_id "
        << "join rule_smiles to_smiles on to_smiles.id = rule.to_smiles_id ";
    if (filter_property) {
        sql << "where property_name.name = ? ";
    }
    sql
        << "order by "
        << "property_name.name, from_smiles.smiles, to_smiles.smiles, "
        << "rule_env.radius, stats.rule_environment_id";
    return sql.str();
}

template <typename Row>
RuleEnvironmentStatistics read_rule_environment_statistics_row(const Row& row) {
    return RuleEnvironmentStatistics(
        row.template GetValue<std::uint64_t>(0),
        row.template GetValue<std::string>(1),
        row.template GetValue<std::string>(2),
        row.template GetValue<std::string>(3),
        static_cast<unsigned int>(row.template GetValue<int>(4)),
        row.template GetValue<std::string>(5),
        row.template GetValue<std::string>(6),
        row.template GetValue<std::string>(7),
        row.template GetValue<std::uint32_t>(8),
        row.template GetValue<double>(9),
        row.template GetValue<bool>(10),
        row.template GetValue<double>(11),
        row.template GetValue<bool>(12),
        row.template GetValue<double>(13),
        row.template GetValue<bool>(14),
        row.template GetValue<double>(15),
        row.template GetValue<double>(16),
        row.template GetValue<double>(17),
        row.template GetValue<double>(18),
        row.template GetValue<double>(19),
        row.template GetValue<double>(20),
        row.template GetValue<bool>(21),
        row.template GetValue<double>(22),
        row.template GetValue<bool>(23),
        row.template GetValue<double>(24)
    );
}

std::vector<RuleEnvironmentStatistics> collect_rule_environment_statistics(
    duckdb::QueryResult& result
) {
    std::vector<RuleEnvironmentStatistics> rows;
    for (const auto& row : result) {
        rows.push_back(read_rule_environment_statistics_row(row));
    }
    return rows;
}

// Read matched pairs from a pair query and attach their numeric properties.
//
// Pairs are materialized first, then every property value for the compounds
// they reference is fetched in a single query and joined in memory. This
// replaces the prior N+1 pattern (one property query per pair) with two
// queries total, which scales far better on large result sets.
std::vector<MatchedPair> read_pairs_with_properties(
    const std::unique_ptr<duckdb::Connection>& connection,
    duckdb::QueryResult& result
) {
    std::vector<MatchedPair> pairs;
    std::set<std::uint64_t> compound_ids;
    for (const auto& row : result) {
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
        compound_ids.insert(pair.GetSourceMoleculeId());
        compound_ids.insert(pair.GetTargetMoleculeId());
        pairs.push_back(std::move(pair));
    }

    if (pairs.empty()) {
        return pairs;
    }

    // Fetch all property values for the referenced compounds in one query,
    // keyed by (compound_id, property_name) for in-memory join below.
    std::ostringstream id_list;
    bool first = true;
    for (const std::uint64_t compound_id : compound_ids) {
        if (!first) {
            id_list << ",";
        }
        id_list << compound_id;
        first = false;
    }

    const std::string property_sql =
        "select compound_property.compound_id, property_name.name, "
        "compound_property.value "
        "from compound_property "
        "join property_name on property_name.id = compound_property.property_name_id "
        "where compound_property.compound_id in (" + id_list.str() + ")";
    std::unique_ptr<duckdb::QueryResult> property_result =
        connection->Query(property_sql);
    if (!property_result) {
        throw StorageError("DuckDB query returned no result: " + property_sql);
    }
    require_success(*property_result, property_sql);

    // property_name -> value for each compound, with names kept in sorted
    // order so SetProperty is invoked in the same order as the prior
    // ``order by property_name.name`` query.
    std::map<std::uint64_t, std::map<std::string, double>> values_by_compound;
    for (const auto& property_row : *property_result) {
        values_by_compound[property_row.GetValue<std::uint64_t>(0)]
            [property_row.GetValue<std::string>(1)] =
            property_row.GetValue<double>(2);
    }

    // Attach a property to a pair only when both compounds carry it, matching
    // the inner-join semantics of the prior per-pair query.
    for (MatchedPair& pair : pairs) {
        const auto source_it = values_by_compound.find(pair.GetSourceMoleculeId());
        const auto target_it = values_by_compound.find(pair.GetTargetMoleculeId());
        if (source_it == values_by_compound.end() ||
            target_it == values_by_compound.end()) {
            continue;
        }
        for (const auto& source_entry : source_it->second) {
            const auto target_value_it = target_it->second.find(source_entry.first);
            if (target_value_it != target_it->second.end()) {
                pair.SetProperty(
                    source_entry.first,
                    source_entry.second,
                    target_value_it->second
                );
            }
        }
    }
    return pairs;
}

// --- Bulk-load per-table Appender helpers. Each streams staged rows through a
// duckdb::Appender in the table's schema column order, then Close() flushes.
// On exception the Appender is destroyed without flushing (RAII), so no
// partial rows leak; the caller owns the surrounding transaction. Values are
// appended as duckdb::Value objects (the Appender's Append(Value) overload) to
// match the codebase's existing typed-value convention.

void AppendMolecules(
    duckdb::Connection& connection,
    const std::vector<MoleculeRecord>& molecules
) {
    if (molecules.empty()) {
        return;
    }
    duckdb::Appender appender(connection, "compound");
    for (const MoleculeRecord& molecule : molecules) {
        appender.BeginRow();
        // compound.id is the verbatim analyzer internal id, not a counter.
        appender.Append(duckdb::Value::UBIGINT(molecule.GetInternalId()));
        appender.Append(string_or_null(molecule.GetExternalId()));
        appender.Append(duckdb::Value(molecule.GetCanonicalSmiles()));  // input_smiles
        appender.Append(duckdb::Value(molecule.GetCanonicalSmiles()));  // clean_smiles
        appender.Append(duckdb::Value::UINTEGER(molecule.GetHeavyAtomCount()));
        appender.EndRow();
    }
    appender.Close();
}

void AppendConstants(
    duckdb::Connection& connection,
    const std::vector<NewConstant>& constants
) {
    if (constants.empty()) {
        return;
    }
    duckdb::Appender appender(connection, "constant_smiles");
    for (const NewConstant& constant : constants) {
        appender.BeginRow();
        appender.Append(duckdb::Value::UBIGINT(constant.id));
        appender.Append(duckdb::Value(constant.smiles));
        appender.EndRow();
    }
    appender.Close();
}

void AppendRuleSmiles(
    duckdb::Connection& connection,
    const std::vector<NewRuleSmiles>& rule_smiles
) {
    if (rule_smiles.empty()) {
        return;
    }
    duckdb::Appender appender(connection, "rule_smiles");
    for (const NewRuleSmiles& row : rule_smiles) {
        appender.BeginRow();
        appender.Append(duckdb::Value::UBIGINT(row.id));
        appender.Append(duckdb::Value(row.smiles));
        // num_heavies is never populated on the current per-row path; keep it
        // NULL for equivalence with get_or_create_named_row_id.
        appender.Append(duckdb::Value(nullptr));
        appender.EndRow();
    }
    appender.Close();
}

void AppendRules(
    duckdb::Connection& connection,
    const std::vector<NewRule>& rules
) {
    if (rules.empty()) {
        return;
    }
    duckdb::Appender appender(connection, "rule");
    for (const NewRule& rule : rules) {
        appender.BeginRow();
        appender.Append(duckdb::Value::UBIGINT(rule.id));
        appender.Append(duckdb::Value::UBIGINT(rule.from_id));
        appender.Append(duckdb::Value::UBIGINT(rule.to_id));
        appender.EndRow();
    }
    appender.Close();
}

void AppendFingerprints(
    duckdb::Connection& connection,
    const std::vector<NewFingerprint>& fingerprints
) {
    if (fingerprints.empty()) {
        return;
    }
    duckdb::Appender appender(connection, "environment_fingerprint");
    for (const NewFingerprint& fingerprint : fingerprints) {
        appender.BeginRow();
        appender.Append(duckdb::Value::UBIGINT(fingerprint.id));
        appender.Append(duckdb::Value(fingerprint.smarts));
        appender.Append(duckdb::Value(fingerprint.pseudo));
        appender.Append(duckdb::Value(fingerprint.parent));
        appender.EndRow();
    }
    appender.Close();
}

void AppendRuleEnvironments(
    duckdb::Connection& connection,
    const std::vector<NewRuleEnvironment>& rule_environments
) {
    if (rule_environments.empty()) {
        return;
    }
    duckdb::Appender appender(connection, "rule_environment");
    for (const NewRuleEnvironment& row : rule_environments) {
        appender.BeginRow();
        appender.Append(duckdb::Value::UBIGINT(row.id));
        appender.Append(duckdb::Value::UBIGINT(row.rule_id));
        appender.Append(duckdb::Value::UBIGINT(row.fingerprint_id));
        appender.Append(duckdb::Value::INTEGER(row.radius));
        appender.Append(duckdb::Value::UINTEGER(
            static_cast<std::uint32_t>(row.num_pairs)));
        appender.EndRow();
    }
    appender.Close();
}

void AppendConstantEnvironments(
    duckdb::Connection& connection,
    const std::vector<NewConstantEnvironment>& rows
) {
    if (rows.empty()) {
        return;
    }
    duckdb::Appender appender(connection, "constant_environment");
    for (const NewConstantEnvironment& row : rows) {
        appender.BeginRow();
        appender.Append(duckdb::Value::UBIGINT(row.constant_id));
        appender.Append(duckdb::Value::INTEGER(row.radius));
        appender.Append(duckdb::Value::UBIGINT(row.fingerprint_id));
        appender.EndRow();
    }
    appender.Close();
}

void AppendPairs(
    duckdb::Connection& connection,
    const std::vector<PairRow>& pair_rows
) {
    if (pair_rows.empty()) {
        return;
    }
    duckdb::Appender appender(connection, "pair");
    for (const PairRow& row : pair_rows) {
        appender.BeginRow();
        appender.Append(duckdb::Value::UBIGINT(row.id));
        appender.Append(duckdb::Value::UBIGINT(row.rule_id));
        appender.Append(duckdb::Value::UBIGINT(row.constant_id));
        appender.Append(duckdb::Value::UBIGINT(row.compound1_id));
        appender.Append(duckdb::Value::UBIGINT(row.compound2_id));
        appender.Append(duckdb::Value::UINTEGER(row.cut_count));
        appender.Append(duckdb::Value::INTEGER(row.heavy_atom_delta));
        appender.Append(duckdb::Value::INTEGER(row.heavy_bond_delta));
        appender.EndRow();
    }
    appender.Close();
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

    // Ensure the schema (and its indexes) exists on open. InitializeSchema is
    // idempotent (every table/index uses ``if not exists``), so this also
    // backfills indexes added in later versions onto previously created
    // database files, not just freshly created ones.
    InitializeSchema();
}

DuckDBStore::~DuckDBStore() = default;

void DuckDBStore::RequireCompatibleSchemaOrThrow() {
    // A versioned store always has a dataset row carrying oemmpa_schema_version;
    // a populated store with a `pair` table but no such row is a pre-versioned
    // legacy store (e.g. one written only via AddPairs, which never inserts the
    // dataset row).
    auto result = connection_->Query(
        "select oemmpa_schema_version from dataset where id = 1");
    require_success(*result, "read schema version");
    for (const auto& row : *result) {
        const std::uint32_t version = row.GetValue<std::uint32_t>(0);
        if (version != kOemmpaSchemaVersion) {
            throw StorageError(
                "DuckDB store was written by oemmpa schema version " +
                std::to_string(version) + " but this build requires version " +
                std::to_string(kOemmpaSchemaVersion) +
                "; rebuild the store from source");
        }
        return;
    }
    throw StorageError(
        "DuckDB store predates schema versioning (no dataset version row); "
        "rebuild the store from source");
}

void DuckDBStore::InitializeSchema() {
    // Gate BEFORE any DDL. An existing store must be validated first, so a
    // legacy store fails with a deterministic rebuild-required error rather than
    // erroring on future version-specific DDL (e.g. columns or indexes added in
    // a later revision).
    const bool fresh_database = !HasTable("pair");
    if (!fresh_database) {
        RequireCompatibleSchemaOrThrow();
    }

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
            "num_rule_environments ubigint,"
            "num_rule_environment_stats ubigint"
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
            "create table if not exists rule_environment_statistics ("
            "id ubigint primary key,"
            "rule_environment_id ubigint not null references rule_environment(id),"
            "property_name_id ubigint not null references property_name(id),"
            "count uinteger not null,"
            "avg double not null,"
            "std double,"
            "kurtosis double,"
            "skewness double,"
            "min double not null,"
            "q1 double not null,"
            "median double not null,"
            "q3 double not null,"
            "max double not null,"
            "paired_t double,"
            "p_value double,"
            "unique (rule_environment_id, property_name_id)"
            ")"
        );
        Execute(
            "create table if not exists constant_smiles ("
            "id ubigint primary key,"
            "smiles varchar not null unique"
            ")"
        );
        Execute(
            "create table if not exists constant_environment ("
            "constant_id ubigint not null references constant_smiles(id),"
            "radius integer not null,"
            "environment_fingerprint_id ubigint not null "
            "references environment_fingerprint(id),"
            "unique (constant_id, radius)"
            ")"
        );
        Execute(
            "create table if not exists pair ("
            "id ubigint primary key,"
            "rule_id ubigint not null references rule(id),"
            "constant_id ubigint not null references constant_smiles(id),"
            "compound1_id ubigint not null references compound(id),"
            "compound2_id ubigint not null references compound(id),"
            "cut_count uinteger not null,"
            "heavy_atom_delta integer not null,"
            "heavy_bond_delta integer not null,"
            "unique (compound1_id, compound2_id, rule_id, constant_id)"
            ")"
        );
        // The pair foreign keys back the hot pair-query joins and the
        // per-pair lookups during bulk loads; the unique constraints on the
        // other tables already provide their backing indexes.
        Execute(
            "create index if not exists pair_rule_idx on pair(rule_id)"
        );
        Execute(
            "create index if not exists pair_constant_idx on pair(constant_id)"
        );
        Execute(
            "create index if not exists pair_compound1_idx on pair(compound1_id)"
        );
        Execute(
            "create index if not exists pair_compound2_idx on pair(compound2_id)"
        );

        // Back each id column with a DuckDB sequence instead of recomputing
        // max(id)+1 on every insert. Each sequence is seeded at max(id)+1 so it
        // is correct on a fresh database (max 0 -> start 1) and on a database
        // whose rows predate the sequence; "if not exists" preserves an
        // already-advanced sequence across reopen.
        for (const char* table : {
            "compound", "compound_property", "constant_smiles",
            "environment_fingerprint", "pair", "property_name",
            "rule", "rule_environment", "rule_smiles"
        }) {
            const std::string table_name = table;
            Execute(
                "create sequence if not exists " + id_sequence_name(table_name) +
                " start " + std::to_string(get_max_id(connection_, table_name) + 1)
            );
        }

        if (fresh_database) {
            // Stamp the version eagerly (not lazily via refresh_dataset_counts) so
            // every store carries its revision even before any counts refresh.
            Execute(
                "insert into dataset (id, oemmpa_schema_version, title, "
                "fragment_options, index_options, is_symmetric) "
                "values (1, " + std::to_string(kOemmpaSchemaVersion) +
                ", '', '', '', true) on conflict (id) do nothing");
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

void DuckDBStore::Execute(const std::string& sql) {
    if (!connection_) {
        throw StorageError("DuckDB connection is not open");
    }
    // A rollback discards rows whose ids may be cached, so invalidate the
    // in-memory id caches to avoid handing out stale ids afterwards.
    if (sql == "rollback") {
        ClearIdCaches();
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

LoadReport DuckDBStore::AddMoleculesFromSmilesFile(
    const std::string& smiles_path,
    const OEDESALT::Desalter* desalter
) {
    InitializeSchema();

    std::ifstream input(smiles_path);
    if (!input) {
        throw StorageError("failed to open SMILES file: " + smiles_path);
    }

    LoadReport report;
    std::uint64_t next_id = get_next_id(connection_, "compound", "id");
    std::string line;
    unsigned int row_number = 0;
    // Wrap the row inserts in a single transaction: per-row failures are
    // absorbed into the report, so a rollback only happens on an unexpected
    // error, and the bulk insert avoids a commit per accepted row.
    Execute("begin transaction");
    try {
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

            MoleculeRecord molecule;
            try {
                molecule = MoleculeRecord::FromSmiles(
                    static_cast<unsigned int>(next_id),
                    smiles,
                    external_id,
                    desalter
                );
                AddMolecule(molecule);
            } catch (const std::exception& exc) {
                report.RecordRejected(row_number, exc.what());
                continue;
            }
            report.RecordAccepted(external_id, molecule.GetStrippedNames());
            ++next_id;
        }
        Execute("commit");
    } catch (...) {
        try {
            Execute("rollback");
        } catch (const StorageError&) {
        }
        throw;
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
    // Insert all property rows in one transaction; per-row failures are
    // recorded in the report rather than aborting the load. The statistics
    // refresh below manages its own transaction, so commit before it runs.
    Execute("begin transaction");
    try {
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
                    throw StorageError(
                        "unknown molecule id in property CSV: " + external_id
                    );
                }

                std::vector<std::pair<std::string, double>> parsed_values;
                for (const std::string& property_name : resolved_property_columns) {
                    const std::string value_text =
                        fields[header_index.at(property_name)];
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
        Execute("commit");
    } catch (...) {
        try {
            Execute("rollback");
        } catch (const StorageError&) {
        }
        throw;
    }

    RefreshRuleEnvironmentStatistics();

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

const std::vector<EnvironmentFingerprint>& DuckDBStore::constant_fingerprints(
    const std::string& constant_smiles
) {
    auto cached = constant_fingerprint_cache_.find(constant_smiles);
    if (cached != constant_fingerprint_cache_.end()) {
        return cached->second;
    }
    auto inserted = constant_fingerprint_cache_.emplace(
        constant_smiles,
        ComputeConstantEnvironmentFingerprints(constant_smiles, 0, 5)
    );
    return inserted.first->second;
}

void DuckDBStore::AddPair(const MatchedPair& pair) {
    // One write implementation: the single-pair overload is the bulk path with a
    // one-element vector. AddPairs owns the transaction and reconciles sequences,
    // so standalone AddPair stays atomic and sequence-consistent.
    AddPairs(std::vector<MatchedPair>{pair});
}

void DuckDBStore::ClearIdCaches() {
    constant_id_cache_.clear();
    rule_smiles_id_cache_.clear();
    rule_id_cache_.clear();
    fingerprint_id_cache_.clear();
    rule_environment_id_cache_.clear();
    pair_identity_cache_.clear();
    constant_environment_ids_.clear();
}

std::uint64_t DuckDBStore::seed_counter(const std::string& table_name) {
    return get_max_id(connection_, table_name) + 1;
}

void DuckDBStore::ReconcileSequences() {
    // DuckDB 1.5.4 does not support ALTER SEQUENCE ... RESTART, so drop and
    // recreate each id sequence starting at max(id)+1. Mirrors the seeding in
    // InitializeSchema so a subsequent legacy nextval insert cannot collide
    // with any bulk- or verbatim-assigned id.
    for (const char* table : {
        "compound", "compound_property", "constant_smiles",
        "environment_fingerprint", "pair", "property_name",
        "rule", "rule_environment", "rule_smiles"
    }) {
        const std::string table_name = table;
        const std::uint64_t start = get_max_id(connection_, table_name) + 1;
        Execute("drop sequence if exists " + id_sequence_name(table_name));
        Execute(
            "create sequence " + id_sequence_name(table_name) +
            " start " + std::to_string(start)
        );
    }
}

std::uint64_t DuckDBStore::cached_named_row_id(
    std::unordered_map<std::string, std::uint64_t>& cache,
    const std::string& table_name,
    const std::string& value
) {
    auto it = cache.find(value);
    if (it != cache.end()) {
        return it->second;
    }
    const std::uint64_t id = get_or_create_named_row_id(
        connection_, table_name, "id", "smiles", value);
    cache.emplace(value, id);
    return id;
}

void DuckDBStore::AddPairs(const std::vector<MatchedPair>& pairs) {
    // Own a transaction only when the caller has not already opened one. This
    // keeps a standalone AddPairs/AddPair atomic (and reconciles id sequences),
    // while still allowing use inside a caller-managed transaction -- DuckDB
    // 1.5.4 errors on a nested "begin transaction", so we must not start one
    // when HasActiveTransaction() is already true. When the caller owns the
    // transaction, it is responsible for commit/rollback; we neither commit nor
    // reconcile sequences here (the owner does that, as SaveTo does).
    const bool owns_transaction = !connection_->HasActiveTransaction();
    if (!owns_transaction) {
        // Molecule existence is enforced by the pair FKs; no molecules here.
        // AppendBulk reconciles id sequences itself, so the caller-owned
        // transaction gets consistent sequences on commit too.
        AppendBulk({}, pairs);
        return;
    }

    Execute("begin transaction");
    try {
        // AppendBulk appends and reconciles id sequences within this transaction.
        AppendBulk({}, pairs);
        Execute("commit");
    } catch (...) {
        try {
            Execute("rollback");
        } catch (const StorageError&) {
        }
        ClearIdCaches();
        throw;
    }
    ClearIdCaches();
}

void DuckDBStore::PreloadIdCaches() {
    // Populate member id caches from existing rows so a non-empty store reuses
    // ids. On a fresh store these queries return nothing (cheap).
    auto load_named = [&](const char* table,
                          std::unordered_map<std::string, std::uint64_t>& cache) {
        auto result = connection_->Query(
            std::string("select smiles, id from ") + table);
        require_success(*result, "preload " + std::string(table));
        for (const auto& row : *result) {
            cache.emplace(row.GetValue<std::string>(0),
                          row.GetValue<std::uint64_t>(1));
        }
    };
    load_named("constant_smiles", constant_id_cache_);
    load_named("rule_smiles", rule_smiles_id_cache_);

    auto rules = connection_->Query(
        "select from_smiles_id, to_smiles_id, id from rule");
    require_success(*rules, "preload rule");
    for (const auto& row : *rules) {
        rule_id_cache_.emplace(
            std::make_pair(row.GetValue<std::uint64_t>(0),
                           row.GetValue<std::uint64_t>(1)),
            row.GetValue<std::uint64_t>(2));
    }

    auto fingerprints = connection_->Query(
        "select smarts, pseudosmiles, parent_smarts, id from environment_fingerprint");
    require_success(*fingerprints, "preload environment_fingerprint");
    for (const auto& row : *fingerprints) {
        const std::string key = row.GetValue<std::string>(0) + "\x1f" +
            row.GetValue<std::string>(1) + "\x1f" + row.GetValue<std::string>(2);
        fingerprint_id_cache_.emplace(key, row.GetValue<std::uint64_t>(3));
    }

    auto environments = connection_->Query(
        "select rule_id, environment_fingerprint_id, radius, id from rule_environment");
    require_success(*environments, "preload rule_environment");
    for (const auto& row : *environments) {
        rule_environment_id_cache_.emplace(
            std::make_tuple(row.GetValue<std::uint64_t>(0),
                            row.GetValue<std::uint64_t>(1),
                            row.GetValue<std::int32_t>(2)),
            row.GetValue<std::uint64_t>(3));
    }

    auto pairs_existing = connection_->Query(
        "select compound1_id, compound2_id, rule_id, constant_id from pair");
    require_success(*pairs_existing, "preload pair identities");
    for (const auto& row : *pairs_existing) {
        pair_identity_cache_.emplace(std::make_tuple(
            row.GetValue<std::uint64_t>(0), row.GetValue<std::uint64_t>(1),
            row.GetValue<std::uint64_t>(2), row.GetValue<std::uint64_t>(3)));
    }
    auto ce_existing = connection_->Query(
        "select distinct constant_id from constant_environment");
    require_success(*ce_existing, "preload constant_environment ids");
    for (const auto& row : *ce_existing) {
        constant_environment_ids_.emplace(row.GetValue<std::uint64_t>(0));
    }
}

void DuckDBStore::AppendBulk(
    const std::vector<MoleculeRecord>& molecules,
    const std::vector<MatchedPair>& pairs
) {
    // --- Seed dimension/pair counters from current max(id). compound is NOT
    // counter-assigned: its id is the analyzer internal id, written verbatim.
    BulkIdCounter constant_counter{seed_counter("constant_smiles")};
    BulkIdCounter rule_smiles_counter{seed_counter("rule_smiles")};
    BulkIdCounter rule_counter{seed_counter("rule")};
    BulkIdCounter fingerprint_counter{seed_counter("environment_fingerprint")};
    BulkIdCounter rule_environment_counter{seed_counter("rule_environment")};
    BulkIdCounter pair_counter{seed_counter("pair")};

    // --- Pre-load existing natural-key -> id rows so a non-empty store reuses
    // ids instead of re-minting them (preserves get_or_create_* semantics).
    // Reuses the member caches; they are cleared by the transaction owner.
    PreloadIdCaches();

    // Validate molecules up front with the same semantics as the legacy
    // per-row AddMolecule, so the bulk path preserves its exception contract:
    // a duplicate external (public) id raises DuplicateIdError -- a distinct,
    // SWIG-exposed type callers catch -- rather than the generic StorageError a
    // raw Appender primary-key violation would surface. Duplicate internal ids
    // and empty canonical SMILES keep their StorageError. Checking before any
    // append means a rejected molecule leaves no partial rows.
    for (const MoleculeRecord& molecule : molecules) {
        if (molecule.GetCanonicalSmiles().empty()) {
            throw StorageError("molecule record has no canonical SMILES");
        }
        if (HasMolecule(molecule.GetInternalId())) {
            throw StorageError(
                "duplicate internal molecule id: " +
                std::to_string(molecule.GetInternalId())
            );
        }
        if (molecule.HasExternalId() &&
            has_external_molecule_id(connection_, molecule.GetExternalId())) {
            throw DuplicateIdError(
                "duplicate molecule id: " + molecule.GetExternalId());
        }
    }

    // Validate every pair's source/target molecule up front, mirroring the
    // legacy per-row AddPair guard, so an orphan-FK pair raises StorageError
    // BEFORE any Appender writes rather than tripping a DuckDB FK constraint
    // mid-append. The latter matters most on a caller-owned transaction: a
    // constraint error aborts the whole active transaction, which would discard
    // the caller's unrelated prior writes. A molecule is valid if it is already
    // persisted OR supplied in this call's molecules batch (the SaveTo case,
    // where compounds are appended in the same transaction just below).
    std::unordered_set<unsigned int> batch_molecule_ids;
    batch_molecule_ids.reserve(molecules.size());
    for (const MoleculeRecord& molecule : molecules) {
        batch_molecule_ids.insert(molecule.GetInternalId());
    }
    const auto require_known_molecule = [&](unsigned int molecule_id,
                                            const char* side) {
        if (batch_molecule_ids.count(molecule_id) == 0 &&
            !HasMolecule(molecule_id)) {
            throw StorageError(
                std::string("cannot add pair with unknown ") + side + " molecule");
        }
    };
    for (const MatchedPair& pair : pairs) {
        require_known_molecule(pair.GetSourceMoleculeId(), "source");
        require_known_molecule(pair.GetTargetMoleculeId(), "target");
    }

    // New-row staging vectors (only rows not already present get appended).
    std::vector<NewConstant> new_constants;
    std::vector<NewRuleSmiles> new_rule_smiles;
    std::vector<NewRule> new_rules;
    std::vector<NewFingerprint> new_fingerprints;
    std::vector<NewRuleEnvironment> new_rule_environments;
    std::vector<NewConstantEnvironment> new_constant_environments;

    // One physical pair row per (compound1, compound2, rule, constant) identity.
    std::vector<PairRow> pair_rows;
    pair_rows.reserve(pairs.size());

    // Local payload map to detect conflicting payloads for the same identity
    // within this batch. The identity cache (pair_identity_cache_) is also seeded
    // from persisted rows, but PreloadIdCaches does not load payloads, so we can
    // only verify payload consistency for identities seen THIS call.
    std::unordered_map<std::tuple<std::uint64_t, std::uint64_t, std::uint64_t,
        std::uint64_t>, std::tuple<unsigned int, int, int>, PairIdentityHash>
        seen_payload;

    // Wrap the whole resolve-then-append body so every failure mode presents the
    // store's StorageError contract to callers. The resolve phase can throw
    // sibling OEMMPAError types that are NOT StorageError -- e.g.
    // EnvironmentFingerprintError from constant_fingerprints() on a malformed
    // constant SMILES -- and the Appender phase can throw a raw duckdb::Exception
    // on a constraint violation. StorageError and DuplicateIdError (the
    // deliberate validation throws below and the duplicate-external-id contract)
    // pass through unchanged; anything else becomes StorageError. The owning
    // transaction (SaveTo / AddPairs) rolls back and clears id caches on throw.
    try {
    for (const MatchedPair& pair : pairs) {
        // Preserve the legacy get_or_create_named_row_id guard: empty normalized
        // values must never be stored. AppendBulk assigns ids directly rather
        // than routing through that helper, so re-check here to keep the bulk
        // path equivalent to legacy AddPair and to avoid persisting empty
        // constant/rule_smiles rows (plus dependent rule/rule_environment/pair
        // rows) for a malformed pair.
        if (pair.GetConstantSmiles().empty()) {
            throw StorageError("cannot store empty normalized value in constant_smiles");
        }
        if (pair.GetSourceVariableSmiles().empty() ||
            pair.GetTargetVariableSmiles().empty()) {
            throw StorageError("cannot store empty normalized value in rule_smiles");
        }
        // constant_smiles
        std::uint64_t constant_id;
        {
            auto it = constant_id_cache_.find(pair.GetConstantSmiles());
            if (it != constant_id_cache_.end()) {
                constant_id = it->second;
            } else {
                constant_id = constant_counter();
                constant_id_cache_.emplace(pair.GetConstantSmiles(), constant_id);
                new_constants.push_back({constant_id, pair.GetConstantSmiles()});
            }
        }
        // rule_smiles (source + target variable)
        auto resolve_rule_smiles = [&](const std::string& smiles) -> std::uint64_t {
            auto it = rule_smiles_id_cache_.find(smiles);
            if (it != rule_smiles_id_cache_.end()) return it->second;
            const std::uint64_t id = rule_smiles_counter();
            rule_smiles_id_cache_.emplace(smiles, id);
            new_rule_smiles.push_back({id, smiles});
            return id;
        };
        const std::uint64_t from_id =
            resolve_rule_smiles(pair.GetSourceVariableSmiles());
        const std::uint64_t to_id =
            resolve_rule_smiles(pair.GetTargetVariableSmiles());
        // rule
        std::uint64_t rule_id;
        {
            const std::pair<std::uint64_t, std::uint64_t> key{from_id, to_id};
            auto it = rule_id_cache_.find(key);
            if (it != rule_id_cache_.end()) {
                rule_id = it->second;
            } else {
                rule_id = rule_counter();
                rule_id_cache_.emplace(key, rule_id);
                new_rules.push_back({rule_id, from_id, to_id});
            }
        }
        // Materialize this constant's per-radius environment memberships once.
        // constant_environment is what the read/statistics joins use to derive
        // each pair's rule_environments from its single physical row, so it must
        // be populated the first time a constant's fingerprints are resolved.
        // Idempotent across reloads via constant_environment_ids_.
        if (constant_environment_ids_.find(constant_id) ==
            constant_environment_ids_.end()) {
            constant_environment_ids_.insert(constant_id);
            for (const EnvironmentFingerprint& fingerprint :
                 constant_fingerprints(pair.GetConstantSmiles())) {
                // fingerprint_id resolution below already caches ids; reuse it.
                const std::string key = fingerprint.GetSmarts() + "\x1f" +
                    fingerprint.GetPseudoSmiles() + "\x1f" +
                    fingerprint.GetParentSmarts();
                std::uint64_t fingerprint_id;
                auto it = fingerprint_id_cache_.find(key);
                if (it != fingerprint_id_cache_.end()) {
                    fingerprint_id = it->second;
                } else {
                    fingerprint_id = fingerprint_counter();
                    fingerprint_id_cache_.emplace(key, fingerprint_id);
                    new_fingerprints.push_back({
                        fingerprint_id, fingerprint.GetSmarts(),
                        fingerprint.GetPseudoSmiles(),
                        fingerprint.GetParentSmarts()});
                }
                new_constant_environments.push_back({
                    constant_id,
                    static_cast<int>(fingerprint.GetRadius()),
                    fingerprint_id});
            }
        }
        // fan out over radius 0..5 fingerprints (cached per constant)
        for (const EnvironmentFingerprint& fingerprint :
             constant_fingerprints(pair.GetConstantSmiles())) {
            std::uint64_t fingerprint_id;
            {
                const std::string key = fingerprint.GetSmarts() + "\x1f" +
                    fingerprint.GetPseudoSmiles() + "\x1f" +
                    fingerprint.GetParentSmarts();
                auto it = fingerprint_id_cache_.find(key);
                if (it != fingerprint_id_cache_.end()) {
                    fingerprint_id = it->second;
                } else {
                    fingerprint_id = fingerprint_counter();
                    fingerprint_id_cache_.emplace(key, fingerprint_id);
                    new_fingerprints.push_back({
                        fingerprint_id, fingerprint.GetSmarts(),
                        fingerprint.GetPseudoSmiles(), fingerprint.GetParentSmarts()});
                }
            }
            // Ensure the (rule, fingerprint, radius) rule_environment row
            // exists. num_pairs is left 0 here and set by the post-append
            // set-based UPDATE, which derives each environment's count from the
            // physical pairs via the constant_environment reconstruction join.
            const std::tuple<std::uint64_t, std::uint64_t, int>
                rule_environment_key{
                    rule_id, fingerprint_id,
                    static_cast<int>(fingerprint.GetRadius())};
            if (rule_environment_id_cache_.find(rule_environment_key) ==
                rule_environment_id_cache_.end()) {
                const std::uint64_t rule_environment_id =
                    rule_environment_counter();
                rule_environment_id_cache_.emplace(
                    rule_environment_key, rule_environment_id);
                new_rule_environments.push_back({
                    rule_environment_id, rule_id, fingerprint_id,
                    static_cast<int>(fingerprint.GetRadius()), 0});
            }
        }
        // One physical pair row per (compound1, compound2, rule, constant)
        // identity. Dedup via the seeded identity cache so a duplicate pair in a
        // later AddPairs call (or a reopened store) is skipped rather than
        // tripping the unique constraint.
        const std::tuple<std::uint64_t, std::uint64_t, std::uint64_t,
            std::uint64_t> identity{
            pair.GetSourceMoleculeId(), pair.GetTargetMoleculeId(),
            rule_id, constant_id};
        const std::tuple<unsigned int, int, int> payload{
            pair.GetCutCount(), pair.GetHeavyAtomDelta(),
            pair.GetHeavyBondDelta()};
        if (pair_identity_cache_.insert(identity).second) {
            // First sight this load AND not persisted: stage the pair and record
            // its payload so later same-identity pairs can be validated.
            pair_rows.push_back({
                pair_counter(), rule_id, constant_id,
                pair.GetSourceMoleculeId(), pair.GetTargetMoleculeId(),
                pair.GetCutCount(), pair.GetHeavyAtomDelta(),
                pair.GetHeavyBondDelta()});
            seen_payload.emplace(identity, payload);
        } else {
            // Identity already known (either staged earlier this call or preloaded
            // from disk). If staged this call, verify payload consistency.
            auto it = seen_payload.find(identity);
            if (it != seen_payload.end() && it->second != payload) {
                // Same identity but DIFFERENT payload within this batch: reject
                // loudly rather than silently dropping the second pair.
                throw StorageError(
                    "conflicting matched-pair payload for identical "
                    "(compound1, compound2, rule, constant): "
                    "cut_count/heavy_atom_delta/heavy_bond_delta mismatch");
            }
            // Otherwise: either exact duplicate (same identity, same payload) or
            // persisted identity we can't verify -> skip silently as before.
        }
    }

    // --- Phase B: append everything in FK-dependency order via Appenders.
    // compound first (verbatim analyzer ids), then dimensions, then pairs.
    AppendMolecules(*connection_, molecules);
    AppendConstants(*connection_, new_constants);
    AppendRuleSmiles(*connection_, new_rule_smiles);
    AppendRules(*connection_, new_rules);
    AppendFingerprints(*connection_, new_fingerprints);
    AppendConstantEnvironments(*connection_, new_constant_environments);
    AppendRuleEnvironments(*connection_, new_rule_environments);
    AppendPairs(*connection_, pair_rows);

    // num_pairs is derived: count physical pairs whose (rule, constant) resolve to
    // each environment's fingerprint at its radius. Set-based, off the per-pair
    // hot path; recompute (not accumulate) so incremental reloads stay correct.
    // MUST run after AppendPairs so newly appended pairs are counted.
    Execute(
        "update rule_environment re set num_pairs = ("
        "select count(*) from pair p "
        "join constant_environment ce "
        "on ce.constant_id = p.constant_id and ce.radius = re.radius "
        "where p.rule_id = re.rule_id "
        "and ce.environment_fingerprint_id = re.environment_fingerprint_id)");

    // Reconcile id sequences to max(id)+1 inside this (possibly caller-owned)
    // transaction, so a later legacy nextval insert cannot collide with the
    // counter-assigned or verbatim ids written above. Doing it here rather than
    // in the transaction owner means EVERY entry point is covered -- the
    // standalone AddPairs/AddPair path, a caller-managed transaction (which
    // cannot reach the private ReconcileSequences itself), and SaveTo. DuckDB
    // 1.5.4 permits the DROP/CREATE SEQUENCE DDL inside an open transaction and
    // the reset persists on commit (verified).
    ReconcileSequences();
    } catch (const StorageError&) {
        throw;
    } catch (const DuplicateIdError&) {
        throw;
    } catch (const std::exception& error) {
        throw StorageError(std::string("DuckDB bulk save failed: ") + error.what());
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

void DuckDBStore::RefreshDatasetCounts() {
    Execute("begin transaction");
    try {
        refresh_dataset_counts(connection_);
        Execute("commit");
    } catch (...) {
        try {
            Execute("rollback");
        } catch (const StorageError&) {
        }
        throw;
    }
}

void DuckDBStore::RefreshRuleEnvironmentStatistics() {
    Execute("begin transaction");
    try {
        refresh_rule_environment_statistics(connection_);
        refresh_dataset_counts(connection_);
        Execute("commit");
    } catch (...) {
        try {
            Execute("rollback");
        } catch (const StorageError&) {
        }
        throw;
    }
}

std::uint64_t DuckDBStore::GetRuleEnvironmentStatisticsCount(
    const std::string& property_name
) const {
    if (property_name.empty()) {
        throw StorageError("rule environment statistics property name must not be empty");
    }

    const std::string sql =
        "select count(*) "
        "from rule_environment_statistics stats "
        "join property_name on property_name.id = stats.property_name_id "
        "where property_name.name = ?";
    duckdb::vector<duckdb::Value> values = {
        duckdb::Value(property_name),
    };
    std::unique_ptr<duckdb::QueryResult> result =
        execute_prepared(connection_, sql, std::move(values));

    for (const auto& row : *result) {
        return row.GetValue<std::uint64_t>(0);
    }
    return 0;
}

std::vector<RuleEnvironmentStatistics> DuckDBStore::GetRuleEnvironmentStatistics() const {
    const std::string sql = build_rule_environment_statistics_query(false);
    std::unique_ptr<duckdb::QueryResult> result = connection_->Query(sql);
    if (!result) {
        throw StorageError("DuckDB query returned no result: " + sql);
    }
    require_success(*result, sql);
    return collect_rule_environment_statistics(*result);
}

std::vector<RuleEnvironmentStatistics> DuckDBStore::GetRuleEnvironmentStatistics(
    const std::string& property_name
) const {
    const std::string sql = build_rule_environment_statistics_query(true);
    std::unique_ptr<duckdb::PreparedStatement> statement = connection_->Prepare(sql);
    if (!statement) {
        throw StorageError("DuckDB prepare returned no statement: " + sql);
    }
    if (statement->HasError()) {
        throw StorageError("DuckDB prepare failed: " + statement->GetError());
    }

    std::unique_ptr<duckdb::QueryResult> result =
        statement->Execute(duckdb::Value(property_name));
    if (!result) {
        throw StorageError("DuckDB query returned no result: " + sql);
    }
    require_success(*result, sql);
    return collect_rule_environment_statistics(*result);
}

DatabaseSummary DuckDBStore::GetSummary(bool recount) const {
    if (recount) {
        return DatabaseSummary(
            GetRowCount("compound"),
            GetRowCount("rule"),
            GetRowCount("pair"),
            GetRowCount("rule_environment"),
            GetRowCount("rule_environment_statistics")
        );
    }

    const std::string sql =
        "select "
        "num_compounds, num_rules, num_pairs, num_rule_environments, "
        "num_rule_environment_stats "
        "from dataset "
        "where id = 1 "
        "and num_compounds is not null "
        "and num_rules is not null "
        "and num_pairs is not null "
        "and num_rule_environments is not null "
        "and num_rule_environment_stats is not null";
    std::unique_ptr<duckdb::QueryResult> result = connection_->Query(sql);
    if (!result) {
        throw StorageError("DuckDB query returned no result: " + sql);
    }
    require_success(*result, sql);

    for (const auto& row : *result) {
        return DatabaseSummary(
            row.GetValue<std::uint64_t>(0),
            row.GetValue<std::uint64_t>(1),
            row.GetValue<std::uint64_t>(2),
            row.GetValue<std::uint64_t>(3),
            row.GetValue<std::uint64_t>(4)
        );
    }

    return GetSummary(true);
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

    std::vector<MatchedPair> pairs =
        read_pairs_with_properties(connection_, *result);
    // Apply the variable-fragment size bounds with the SAME predicate as the
    // in-memory path so the store-read and Analyzer query results match. Done
    // before scoring, which collapses candidate groups.
    pairs = filter_pairs_by_variable_bounds(connection_, std::move(pairs), options);
    return apply_pair_scoring(pairs, options.GetScoringOptions());
}

std::vector<MatchedPair> DuckDBStore::GetPairsForRuleEnvironment(
    std::uint64_t rule_environment_id
) const {
    if (rule_environment_id == 0) {
        throw StorageError("rule environment id must be greater than zero");
    }

    const QueryOptions options;
    const std::string sql = build_pair_query(options, rule_environment_id);
    std::unique_ptr<duckdb::QueryResult> result = connection_->Query(sql);
    if (!result) {
        throw StorageError("DuckDB query returned no result: " + sql);
    }
    require_success(*result, sql);

    return read_pairs_with_properties(connection_, *result);
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
