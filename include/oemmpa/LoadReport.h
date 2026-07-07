#ifndef OEMMPA_LOAD_REPORT_H
#define OEMMPA_LOAD_REPORT_H

#include <string>
#include <vector>

namespace OEMMPA {

struct LoadError {
    unsigned int row = 0;
    std::string message;
};

struct AcceptedMolecule {
    std::string external_id;
    std::vector<std::string> stripped_names;
};

class LoadReport {
public:
    void RecordAccepted(const std::string& external_id, const std::vector<std::string>& stripped_names = {});
    void RecordRejected(unsigned int row, const std::string& message);

    unsigned int GetAcceptedCount() const;
    unsigned int GetRejectedCount() const;
    const std::vector<std::string>& GetAcceptedIds() const;
    const std::vector<AcceptedMolecule>& GetAcceptedMolecules() const;
    const std::vector<LoadError>& GetErrors() const;

private:
    std::vector<std::string> accepted_ids_;
    std::vector<AcceptedMolecule> accepted_molecules_;
    std::vector<LoadError> errors_;
};

}  // namespace OEMMPA

#endif  // OEMMPA_LOAD_REPORT_H
