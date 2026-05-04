#ifndef OEMMPA_QUERY_ENVIRONMENT_H
#define OEMMPA_QUERY_ENVIRONMENT_H

#include <string>
#include <vector>

namespace OEMMPA {

/// \brief Local environment observed while fragmenting a query molecule.
class QueryEnvironment {
public:
    QueryEnvironment() = default;
    QueryEnvironment(
        std::string constant_smiles,
        std::string variable_smiles,
        unsigned int cut_count,
        unsigned int radius,
        std::string smarts,
        std::string pseudo_smiles,
        std::string parent_smarts
    );

    const std::string& GetConstantSmiles() const;
    const std::string& GetVariableSmiles() const;
    unsigned int GetCutCount() const;
    unsigned int GetRadius() const;
    const std::string& GetSmarts() const;
    const std::string& GetPseudoSmiles() const;
    const std::string& GetParentSmarts() const;

private:
    std::string constant_smiles_;
    std::string variable_smiles_;
    unsigned int cut_count_ = 0;
    unsigned int radius_ = 0;
    std::string smarts_;
    std::string pseudo_smiles_;
    std::string parent_smarts_;
};

std::vector<QueryEnvironment> ComputeQueryEnvironments(
    const std::string& smiles,
    unsigned int min_radius = 0,
    unsigned int max_radius = 5
);

bool SmilesContainsSubstructure(
    const std::string& smiles,
    const std::string& smarts
);

}  // namespace OEMMPA

#endif  // OEMMPA_QUERY_ENVIRONMENT_H
