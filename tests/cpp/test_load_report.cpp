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

}  // namespace test
}  // namespace OEMMPA
