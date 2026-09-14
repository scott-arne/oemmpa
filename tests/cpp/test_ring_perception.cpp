#include <gtest/gtest.h>

#include "oemmpa/Error.h"
#include "oemmpa/FragmentSmiles.h"
#include "oemmpa/RingPerception.h"
#include "ring_search_bound_guard.h"

#include <oechem.h>

#include <string>
#include <type_traits>

namespace OEMMPA {
namespace test {
namespace {

OEChem::OEGraphMol MolFromSmiles(const std::string& smiles) {
    OEChem::OEGraphMol mol;
    EXPECT_TRUE(OEChem::OESmilesToMol(mol, smiles)) << smiles;
    return mol;
}

/// Decalin: two fused six-rings, so the search runs a breadth-first ring hunt
/// and a forced bound of 1 aborts it on the first queued atom.
const char* const DECALIN = "C1CCC2CCCCC2C1";

/// N-methylpiperidine written with a tagged nitrogen. Rendering consults ring
/// perception only when the cleanup has a tagged nitrogen to judge
/// (build_stereo_cleanup_context asks find_bridgeheads for nothing else), so
/// this is the smallest render that goes through the ring search.
const char* const TAGGED_RING_NITROGEN = "C[N@]1CCCCC1";

}  // namespace

// The port used to answer a failed search with upstream's approximate ring
// finder, which silently changed the bridgehead set the stereo cleanup reads.
// RDKit throws here, its sanitization fails, and Chem.MolFromSmiles returns
// None -- so the error type is the one the fragmenter already treats as "this
// molecule produces nothing".
TEST(RingPerception, RingPerceptionErrorIsAFragmentationError) {
    static_assert(
        std::is_base_of<FragmentationError, RingPerceptionError>::value,
        "a structure the ring search gives up on is a fragmentation failure");

    try {
        throw RingPerceptionError("ring perception gave up");
    } catch (const FragmentationError& error) {
        EXPECT_STREQ(error.what(), "ring perception gave up");
        return;
    }
    FAIL() << "RingPerceptionError did not derive from FragmentationError";
}

TEST(RingPerception, ExceedingTheSearchBoundThrowsInsteadOfApproximating) {
    const OEChem::OEGraphMol mol = MolFromSmiles(DECALIN);

    const ScopedRingSearchQueueBound forced(1);
    EXPECT_THROW(symmetrized_sssr(mol), RingPerceptionError);
    EXPECT_THROW(bridgehead_atoms(mol), RingPerceptionError);
}

// The message has to name the bound, because the only way to act on it is to
// know that RDKit would have rejected the structure outright.
TEST(RingPerception, TheAbortMessageNamesTheBoundAndTheConsequence) {
    const OEChem::OEGraphMol mol = MolFromSmiles(DECALIN);

    const ScopedRingSearchQueueBound forced(1);
    try {
        symmetrized_sssr(mol);
        FAIL() << "the forced bound did not abort the search";
    } catch (const RingPerceptionError& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("MAX_BFSQ_SIZE"), std::string::npos) << message;
        EXPECT_NE(message.find("RDKit"), std::string::npos) << message;
    }
}

TEST(RingPerception, RenderingATaggedRingNitrogenPropagatesTheAbort) {
    OEChem::OEGraphMol mol = MolFromSmiles(TAGGED_RING_NITROGEN);

    const ScopedRingSearchQueueBound forced(1);
    EXPECT_THROW(render_fragment_smiles(mol), RingPerceptionError);
}

// The control for all three: nothing above is a property of these molecules.
// Under the shipped bound the same calls answer normally, which is also what
// makes the forced-bound tests evidence about the bound.
TEST(RingPerception, TheSameMoleculesPerceiveNormallyUnderTheShippedBound) {
    EXPECT_EQ(detail::ring_search_queue_bound(), 200000u);

    const OEChem::OEGraphMol decalin = MolFromSmiles(DECALIN);
    EXPECT_EQ(symmetrized_sssr(decalin).size(), 2u);
    EXPECT_TRUE(bridgehead_atoms(decalin).empty());

    OEChem::OEGraphMol tagged = MolFromSmiles(TAGGED_RING_NITROGEN);
    EXPECT_FALSE(render_fragment_smiles(tagged).empty());
}

// The seam is process-wide, so the guard's restore is load-bearing for every
// other test in this binary.
TEST(RingPerception, TheGuardRestoresTheBoundOnScopeExit) {
    const std::size_t shipped = detail::ring_search_queue_bound();
    {
        const ScopedRingSearchQueueBound forced(1);
        EXPECT_EQ(detail::ring_search_queue_bound(), 1u);
    }
    EXPECT_EQ(detail::ring_search_queue_bound(), shipped);
}

// perceives_rings() is the strict policy's whole test, so it has to answer
// exactly what symmetrized_sssr() would do: true where it returns, false where
// it throws, and nothing in between. An acyclic molecule never starts the ring
// search, so the forced bound cannot touch it -- that is the case that keeps
// the strict policy from rejecting every molecule when a test forces the bound.
TEST(RingPerception, PerceivesRingsAnswersWhatSymmetrizedSssrWouldDo) {
    OEChem::OEGraphMol benzene;
    ASSERT_TRUE(OEChem::OESmilesToMol(benzene, "c1ccccc1"));
    OEChem::OEGraphMol ethanol;
    ASSERT_TRUE(OEChem::OESmilesToMol(ethanol, "CCO"));

    EXPECT_TRUE(perceives_rings(benzene));
    EXPECT_TRUE(perceives_rings(ethanol));

    const ScopedRingSearchQueueBound forced(1);
    EXPECT_THROW(symmetrized_sssr(benzene), RingPerceptionError);
    EXPECT_FALSE(perceives_rings(benzene));
    EXPECT_TRUE(perceives_rings(ethanol));
}

}  // namespace test
}  // namespace OEMMPA
