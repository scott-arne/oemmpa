#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "oemmpa/SeriesNetwork.h"
#include "oemmpa/SeriesQuery.h"
#include "series_network_support.h"

using namespace OEMMPA;
using OEMMPA::test::EdgeKey;
using OEMMPA::test::EdgeKeys;
using OEMMPA::test::ExpectSameNetwork;
using OEMMPA::test::FindEdge;

namespace {

// R-group keys of the hand-built profiles; the loader ranks real ones by SMILES.
const std::uint32_t METHYL = 1;
const std::uint32_t AMINO = 2;
const std::uint32_t HYDROXY = 3;
const std::uint32_t FLUORO = 4;
const std::uint32_t ETHANO = 5;
const std::uint32_t PROPANO = 6;
const std::uint32_t BUTANO = 7;

SeriesNetworkProfile Profile(std::uint32_t cut_count, std::uint64_t id,
                             const std::vector<std::pair<std::uint32_t, double>>& members) {
    SeriesNetworkProfile profile;
    profile.cut_count = cut_count;
    profile.id = id;
    profile.label = "profile " + std::to_string(id);
    for (const auto& [rgroup, value] : members) {
        profile.members.push_back({rgroup, 1U, value});
    }
    return profile;
}

// Seven profiles whose every pairwise measure the tests work out: S1 and S2 value
// methyl, amino, hydroxy and fluoro (1, 2, 3, 4); S3 reverses them (4, 3, 2, 1);
// S4 is S1 without fluoro; S5 is flat at 1; S6 and S7 are two-cut profiles valued
// (1, 2, 3) over three linkers. Pairwise, over the shared R-groups:
//   S1-S2 rho 1, crmsd 0 (overlap 4)      S1-S3 rho -1, crmsd sqrt 5 (4)
//   S1-S4 rho 1, crmsd 0 (3)              S1-S5 undefined, crmsd sqrt 1.25 (4)
//   S2 pairs as S1's                      S3-S4 rho -1, crmsd sqrt(8/3) (3)
//   S3-S5 undefined, crmsd sqrt 1.25 (4)  S4-S5 undefined, crmsd sqrt(2/3) (3)
//   S6-S7 rho 1, crmsd 0 (3)
std::vector<SeriesNetworkProfile> HandBuiltProfiles() {
    return {
        Profile(1, 1, {{METHYL, 1.0}, {AMINO, 2.0}, {HYDROXY, 3.0}, {FLUORO, 4.0}}),
        Profile(1, 2, {{METHYL, 1.0}, {AMINO, 2.0}, {HYDROXY, 3.0}, {FLUORO, 4.0}}),
        Profile(1, 3, {{METHYL, 4.0}, {AMINO, 3.0}, {HYDROXY, 2.0}, {FLUORO, 1.0}}),
        Profile(1, 4, {{METHYL, 1.0}, {AMINO, 2.0}, {HYDROXY, 3.0}}),
        Profile(1, 5, {{METHYL, 1.0}, {AMINO, 1.0}, {HYDROXY, 1.0}, {FLUORO, 1.0}}),
        Profile(2, 6, {{ETHANO, 1.0}, {PROPANO, 2.0}, {BUTANO, 3.0}}),
        Profile(2, 7, {{ETHANO, 1.0}, {PROPANO, 2.0}, {BUTANO, 3.0}}),
    };
}

SeriesNetworkRequest Request(SimilarityMeasure measure, std::uint32_t top) {
    SeriesNetworkRequest request;
    request.property_name = "act";
    request.measure = measure;
    request.top = top;
    return request;
}

}  // namespace

TEST(SeriesNetworkScoring, SpearmanTopOneKeepsEachSourcesBestAndUnionsThem) {
    const SeriesNetworkResult result =
        score_series_network(HandBuiltProfiles(), Request(SimilarityMeasure::SPEARMAN, 1));
    // S1 and S2 pick each other (rho 1 over four beats S4's rho 1 over three); S3
    // picks S1 (rho -1 ties S2 and S4; overlap 4 beats S4's 3, the lower id beats
    // S2); S4 picks S1 (rho 1 ties S2; lower id); S5 has only undefined
    // correlations and picks S1 (overlap 4, lowest id); S6 and S7 pick each other.
    // S1 ends with four incident edges although top is 1: the others chose it.
    EXPECT_EQ(EdgeKeys(result), (std::vector<EdgeKey>{{1, 1, 2}, {1, 1, 3}, {1, 1, 4}, {1, 1, 5}, {2, 6, 7}}));
    const SeriesNetworkEdge* reverse = FindEdge(result, 1, 1, 3);
    ASSERT_NE(reverse, nullptr);
    EXPECT_EQ(reverse->overlap, 4U);
    EXPECT_TRUE(reverse->has_spearman);
    EXPECT_DOUBLE_EQ(reverse->spearman, -1.0);
    EXPECT_NEAR(reverse->centred_rmsd, std::sqrt(5.0), 1e-12);
    const SeriesNetworkEdge* flat = FindEdge(result, 1, 1, 5);
    ASSERT_NE(flat, nullptr);
    EXPECT_FALSE(flat->has_spearman);
    EXPECT_NEAR(flat->centred_rmsd, std::sqrt(1.25), 1e-12);
    EXPECT_EQ(result.environment_radius, -1);
    // Every profile is a node, by (cut_count, id), sized by its valued R-groups.
    ASSERT_EQ(result.nodes.size(), 7U);
    EXPECT_EQ(result.nodes[3].id, 4U);
    EXPECT_EQ(result.nodes[3].size, 3U);
    EXPECT_EQ(result.nodes[3].label, "profile 4");
    EXPECT_EQ(result.nodes[5].cut_count, 2U);
}

TEST(SeriesNetworkScoring, CentredRmsdTopOneOrdersAscending) {
    const SeriesNetworkResult result =
        score_series_network(HandBuiltProfiles(), Request(SimilarityMeasure::CENTRED_RMSD, 1));
    // S1 and S2 are 0 apart; S3 is nearest S5 (sqrt 1.25 against sqrt 5 and S4's
    // sqrt(8/3)); S4 picks S1 (0, ties S2, lower id); S5 is nearest S4 (sqrt(2/3)).
    EXPECT_EQ(EdgeKeys(result), (std::vector<EdgeKey>{{1, 1, 2}, {1, 1, 4}, {1, 3, 5}, {1, 4, 5}, {2, 6, 7}}));
}

TEST(SeriesNetworkScoring, SpearmanThresholdKeepsDefinedCorrelationsAtOrAboveIt) {
    SeriesNetworkRequest request = Request(SimilarityMeasure::SPEARMAN, 10);
    request.has_threshold = true;
    request.threshold = 0.5;
    const SeriesNetworkResult result = score_series_network(HandBuiltProfiles(), request);
    // Only the rho 1 pairs pass: S3's -1 does not, S5's undefined correlations never do.
    EXPECT_EQ(EdgeKeys(result), (std::vector<EdgeKey>{{1, 1, 2}, {1, 1, 4}, {1, 2, 4}, {2, 6, 7}}));
    request.threshold = 1.0;   // Inclusive.
    EXPECT_EQ(EdgeKeys(score_series_network(HandBuiltProfiles(), request)), EdgeKeys(result));
}

TEST(SeriesNetworkScoring, CentredRmsdThresholdKeepsPairsAtOrBelowIt) {
    SeriesNetworkRequest request = Request(SimilarityMeasure::CENTRED_RMSD, 10);
    request.has_threshold = true;
    request.threshold = 1.0;
    EXPECT_EQ(EdgeKeys(score_series_network(HandBuiltProfiles(), request)),
              (std::vector<EdgeKey>{{1, 1, 2}, {1, 1, 4}, {1, 2, 4}, {1, 4, 5}, {2, 6, 7}}));
    request.threshold = 0.0;   // Inclusive: only identical centred profiles.
    EXPECT_EQ(EdgeKeys(score_series_network(HandBuiltProfiles(), request)),
              (std::vector<EdgeKey>{{1, 1, 2}, {1, 1, 4}, {1, 2, 4}, {2, 6, 7}}));
}

TEST(SeriesNetworkScoring, MemberBoundAndMinimumLengthDecideEligibility) {
    std::vector<SeriesNetworkProfile> profiles = HandBuiltProfiles();
    for (SeriesNetworkProfile& profile : profiles) {
        for (SeriesNetworkMember& member : profile.members) {
            if (member.rgroup == FLUORO) {
                member.num_heavies = 5;   // Fluoro stands in for a heavy R-group.
            }
        }
    }
    SeriesNetworkRequest request = Request(SimilarityMeasure::SPEARMAN, 10);
    request.selection.max_variable_heavy_atoms = 4;
    const SeriesNetworkResult bounded = score_series_network(profiles, request);
    // Without fluoro every profile holds three valued R-groups.
    ASSERT_EQ(bounded.nodes.size(), 7U);
    for (const SeriesNetworkNode& node : bounded.nodes) {
        EXPECT_EQ(node.size, 3U) << node.id;
    }
    const SeriesNetworkEdge* reverse = FindEdge(bounded, 1, 1, 3);
    ASSERT_NE(reverse, nullptr);
    EXPECT_EQ(reverse->overlap, 3U);
    // min_series_length 4 then leaves no profile at all.
    request.selection.min_series_length = 4;
    const SeriesNetworkResult strict = score_series_network(profiles, request);
    EXPECT_TRUE(strict.nodes.empty());
    EXPECT_TRUE(strict.edges.empty());
    // A profile with fewer than two valued R-groups is never a node, whatever
    // min_series_length says.
    request.selection = SeriesSelection();
    request.selection.min_series_length = 0;
    const std::vector<SeriesNetworkProfile> single = {
        Profile(1, 9, {{METHYL, 1.0}}), Profile(1, 10, {{METHYL, 2.0}, {AMINO, 1.0}})};
    const SeriesNetworkResult lone = score_series_network(single, request);
    ASSERT_EQ(lone.nodes.size(), 1U);
    EXPECT_EQ(lone.nodes[0].id, 10U);
}

TEST(SeriesNetworkScoring, MinOverlapCountsSharedValuedRGroups) {
    SeriesNetworkRequest request = Request(SimilarityMeasure::SPEARMAN, 10);
    request.min_overlap = 4;
    // Only the four-member profiles S1, S2, S3 and S5 share four R-groups.
    EXPECT_EQ(EdgeKeys(score_series_network(HandBuiltProfiles(), request)),
              (std::vector<EdgeKey>{{1, 1, 2}, {1, 1, 3}, {1, 1, 5}, {1, 2, 3}, {1, 2, 5}, {1, 3, 5}}));
}

TEST(SeriesNetworkScoring, ResultsAreIdenticalForOneTwoAndEightThreads) {
    // Forty profiles over twelve R-groups, valued by a fixed rule so ties occur.
    std::vector<SeriesNetworkProfile> profiles;
    for (std::uint64_t id = 1; id <= 40; ++id) {
        std::vector<std::pair<std::uint32_t, double>> members;
        for (std::uint32_t rgroup = 1; rgroup <= 12; ++rgroup) {
            if ((id * 7 + rgroup * 3) % 4 != 0) {
                members.emplace_back(rgroup, static_cast<double>((id * rgroup * 13) % 17) / 4.0);
            }
        }
        profiles.push_back(Profile(1, id, members));
    }
    for (const SimilarityMeasure measure : {SimilarityMeasure::SPEARMAN, SimilarityMeasure::CENTRED_RMSD}) {
        SeriesNetworkRequest request = Request(measure, 3);
        request.threads = 1;
        const SeriesNetworkResult serial = score_series_network(profiles, request);
        ASSERT_FALSE(serial.edges.empty());
        for (const std::uint32_t threads : {2U, 8U}) {
            request.threads = threads;
            ExpectSameNetwork(serial, score_series_network(profiles, request));
        }
    }
}

TEST(SeriesNetworkScoring, MeasuresAreSymmetricBitForBit) {
    // The union keeps one copy of a pair both endpoints selected; that is only a
    // union if both directions give the same values to the last bit.
    const std::vector<std::vector<double>> vectors = {
        {1.0, 2.5, 2.5, 7.0}, {0.3, 0.1, 0.2, 0.4}, {3.0, 3.0, 1.0, 2.0}, {1e9, 1.0, -4.0, 0.5}};
    for (const std::vector<double>& x : vectors) {
        for (const std::vector<double>& y : vectors) {
            const RankCorrelation forward = spearman_correlation(x, y);
            const RankCorrelation backward = spearman_correlation(y, x);
            EXPECT_EQ(forward.defined, backward.defined);
            EXPECT_EQ(forward.rho, backward.rho);
            EXPECT_EQ(centred_rmsd(x, y), centred_rmsd(y, x));
        }
    }
}
