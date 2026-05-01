#ifndef OEMMPA_ERROR_H
#define OEMMPA_ERROR_H

#include <stdexcept>
#include <string>

namespace OEMMPA {

class OEMMPAError : public std::runtime_error {
public:
    explicit OEMMPAError(const std::string& message)
        : std::runtime_error(message) {}
};

class InvalidMoleculeError : public OEMMPAError {
public:
    explicit InvalidMoleculeError(const std::string& message)
        : OEMMPAError(message) {}
};

class DuplicateIdError : public OEMMPAError {
public:
    explicit DuplicateIdError(const std::string& message)
        : OEMMPAError(message) {}
};

class MissingPropertyError : public OEMMPAError {
public:
    explicit MissingPropertyError(const std::string& message)
        : OEMMPAError(message) {}
};

class FragmentationError : public OEMMPAError {
public:
    explicit FragmentationError(const std::string& message)
        : OEMMPAError(message) {}
};

class InvalidQueryError : public OEMMPAError {
public:
    explicit InvalidQueryError(const std::string& message)
        : OEMMPAError(message) {}
};

class AnalysisStateError : public OEMMPAError {
public:
    explicit AnalysisStateError(const std::string& message)
        : OEMMPAError(message) {}
};

class StorageError : public OEMMPAError {
public:
    explicit StorageError(const std::string& message)
        : OEMMPAError(message) {}
};

}  // namespace OEMMPA

#endif  // OEMMPA_ERROR_H
