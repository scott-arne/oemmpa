#include <gtest/gtest.h>

#include "oemmpa/Error.h"

namespace OEMMPA {
namespace test {

TEST(ErrorTest, DomainErrorsCarryMessages) {
    InvalidMoleculeError error("invalid molecule: row 3");
    EXPECT_STREQ(error.what(), "invalid molecule: row 3");
}

TEST(ErrorTest, SpecificErrorsDeriveFromOEMMPAError) {
    try {
        throw DuplicateIdError("duplicate external id: cmpd-1");
    } catch (const OEMMPAError& error) {
        EXPECT_STREQ(error.what(), "duplicate external id: cmpd-1");
        return;
    }
    FAIL() << "DuplicateIdError did not derive from OEMMPAError";
}

}  // namespace test
}  // namespace OEMMPA
