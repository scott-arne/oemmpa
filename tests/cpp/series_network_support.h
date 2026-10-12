#pragma once

/// \file
/// Test-only helpers shared by the network tests of test_series_network.cpp and
/// test_series_query.cpp: edge keys and a field-by-field comparison of two
/// network results (the records have no operator==).

#include <gtest/gtest.h>

#include "oemmpa/SeriesQuery.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <vector>

namespace OEMMPA {
namespace test {

/// \brief (cut_count, source_id, target_id) of an undirected edge.
using EdgeKey = std::tuple<std::uint32_t, std::uint64_t, std::uint64_t>;

/// \brief The key of the edge between two profiles of one cut count, lower id first.
inline EdgeKey Key(std::uint32_t cut_count, std::uint64_t a, std::uint64_t b) {
    return EdgeKey{cut_count, std::min(a, b), std::max(a, b)};
}

/// \brief Two doubles are the same result: equal, or both NaN (an overflowed centred RMSD).
inline bool SameDouble(double a, double b) {
    return a == b || (std::isnan(a) && std::isnan(b));
}

/// \brief The keys of a result's edges, in the result's order.
inline std::vector<EdgeKey> EdgeKeys(const SeriesNetworkResult& result) {
    std::vector<EdgeKey> keys;
    for (const SeriesNetworkEdge& edge : result.edges) {
        keys.emplace_back(edge.cut_count, edge.source_id, edge.target_id);
    }
    return keys;
}

/// \brief The edge between two profiles of one cut count, or nullptr.
inline const SeriesNetworkEdge* FindEdge(const SeriesNetworkResult& result, std::uint32_t cut_count,
                                         std::uint64_t a, std::uint64_t b) {
    const EdgeKey wanted = Key(cut_count, a, b);
    for (const SeriesNetworkEdge& edge : result.edges) {
        if (EdgeKey{edge.cut_count, edge.source_id, edge.target_id} == wanted) {
            return &edge;
        }
    }
    return nullptr;
}

/// \brief Expect two results to agree on every field, doubles exactly.
inline void ExpectSameNetwork(const SeriesNetworkResult& expected, const SeriesNetworkResult& actual) {
    EXPECT_EQ(actual.environment_radius, expected.environment_radius);
    ASSERT_EQ(actual.nodes.size(), expected.nodes.size());
    for (std::size_t i = 0; i < expected.nodes.size(); ++i) {
        EXPECT_EQ(actual.nodes[i].id, expected.nodes[i].id) << "node " << i;
        EXPECT_EQ(actual.nodes[i].label, expected.nodes[i].label) << "node " << i;
        EXPECT_EQ(actual.nodes[i].cut_count, expected.nodes[i].cut_count) << "node " << i;
        EXPECT_EQ(actual.nodes[i].size, expected.nodes[i].size) << "node " << i;
    }
    ASSERT_EQ(actual.edges.size(), expected.edges.size());
    for (std::size_t i = 0; i < expected.edges.size(); ++i) {
        const SeriesNetworkEdge& a = actual.edges[i];
        const SeriesNetworkEdge& e = expected.edges[i];
        EXPECT_EQ((EdgeKey{a.cut_count, a.source_id, a.target_id}), (EdgeKey{e.cut_count, e.source_id, e.target_id}))
            << "edge " << i;
        EXPECT_EQ(a.overlap, e.overlap) << "edge " << i;
        EXPECT_EQ(a.has_spearman, e.has_spearman) << "edge " << i;
        EXPECT_TRUE(SameDouble(a.spearman, e.spearman)) << "edge " << i;
        EXPECT_TRUE(SameDouble(a.centred_rmsd, e.centred_rmsd)) << "edge " << i;
    }
}

}  // namespace test
}  // namespace OEMMPA
