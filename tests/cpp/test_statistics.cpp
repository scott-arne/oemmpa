#include <gtest/gtest.h>

#include "oemmpa/Statistics.h"

TEST(Statistics, IncompleteBetaMatchesKnownValues) {
    EXPECT_NEAR(OEMMPA::regularized_incomplete_beta(0.5, 0.5, 0.0), 0.0, 1e-12);
    EXPECT_NEAR(OEMMPA::regularized_incomplete_beta(0.5, 0.5, 1.0), 1.0, 1e-12);
    // I(0.5; 1, 1) = 0.5 exactly.
    EXPECT_NEAR(OEMMPA::regularized_incomplete_beta(1.0, 1.0, 0.5), 0.5, 1e-12);
    // Symmetry: I(x; a, b) = 1 - I(1-x; b, a).
    EXPECT_NEAR(
        OEMMPA::regularized_incomplete_beta(2.0, 3.0, 0.4),
        1.0 - OEMMPA::regularized_incomplete_beta(3.0, 2.0, 0.6),
        1e-12);
}

TEST(Statistics, TwoSidedPValueIsSymmetricAndBounded) {
    EXPECT_NEAR(OEMMPA::two_sided_t_p_value(0.0, 10.0), 1.0, 1e-12);
    EXPECT_NEAR(
        OEMMPA::two_sided_t_p_value(2.5, 10.0),
        OEMMPA::two_sided_t_p_value(-2.5, 10.0),
        1e-15);
    const double p = OEMMPA::two_sided_t_p_value(2.5, 10.0);
    EXPECT_GT(p, 0.0);
    EXPECT_LT(p, 1.0);
}

TEST(Statistics, TwoSidedPValueMatchesScipyReferenceValues) {
    // scipy.stats.t.sf(abs(t), df) * 2, computed once and vendored.
    EXPECT_NEAR(OEMMPA::two_sided_t_p_value(2.228138851986273, 10.0), 0.05, 1e-9);
    EXPECT_NEAR(OEMMPA::two_sided_t_p_value(1.9599639845400545, 1.0e6), 0.05, 1e-5);
}
