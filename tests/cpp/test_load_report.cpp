#include <gtest/gtest.h>

#include "oemmpa/LoadReport.h"

namespace OEMMPA {
namespace test {

TEST(LoadReportTest, TracksAcceptedAndRejectedRows) {
    LoadReport report;
    report.RecordAccepted("cmpd-1");
    report.RecordRejected(4, "invalid molecule");

    EXPECT_EQ(report.GetAcceptedCount(), 1);
    EXPECT_EQ(report.GetRejectedCount(), 1);
    ASSERT_EQ(report.GetAcceptedIds().size(), 1);
    EXPECT_EQ(report.GetAcceptedIds()[0], "cmpd-1");
    ASSERT_EQ(report.GetErrors().size(), 1);
    EXPECT_EQ(report.GetErrors()[0].row, 4);
    EXPECT_EQ(report.GetErrors()[0].message, "invalid molecule");
}

TEST(LoadReport, RecordsStrippedNamesPerAcceptedRow) {
    LoadReport report;
    report.RecordAccepted("mol1", {"Halides"});
    report.RecordAccepted("mol2", {});
    ASSERT_EQ(report.GetAcceptedMolecules().size(), 2u);
    EXPECT_EQ(report.GetAcceptedMolecules()[0].external_id, "mol1");
    EXPECT_EQ(report.GetAcceptedMolecules()[0].stripped_names.size(), 1u);
    EXPECT_TRUE(report.GetAcceptedMolecules()[1].stripped_names.empty());
    // Backward-compatible id accessor still works.
    ASSERT_EQ(report.GetAcceptedIds().size(), 2u);
    EXPECT_EQ(report.GetAcceptedIds()[0], "mol1");
}

}  // namespace test
}  // namespace OEMMPA
