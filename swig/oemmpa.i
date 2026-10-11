// swig/oemmpa.i
// SWIG interface file for oemmpa Python bindings
// The docstring is the reading rule for everything -doxygen copies out of the
// headers: a :raises: line below names the C++ type the header declares, which
// is never what Python catches. The generated module is rebuilt from this file
// on every configure, so the note cannot live in the generated module itself.
%module(docstring="Low-level OEMMPA bindings, generated from the C++ headers. The :raises: lines in these docstrings are copied from those headers and name C++ exception types, such as RingPerceptionError or InvalidQueryError; the classes of those names in this module are wrapped C++ classes, not Python exceptions, and cannot be caught. Every C++ exception raised through these bindings reaches Python as RuntimeError. The handwritten functions of the oemmpa package are the exception: they convert the cases their own docstrings document into ValueError.", "threads"=1) _oemmpa

// ============================================================================
// GIL release configuration (MUST be here, before any wrapper-generating directives)
// ============================================================================
// Release the GIL only around the long, Python-free C++ calls; keep the many
// tiny getters GIL-held. In SWIG 4.4 the executable feature is `nothreadallow`
// (%threadallow == %feature("nothreadallow","0")).
// This global default turns release OFF for ALL wrappers.
%feature("nothreadallow", "1");

%{
#include "oemmpa/oemmpa.h"
#include "oedesalt/Desalter.h"
#include <oechem.h>

using namespace OEMMPA;
%}

// ============================================================================
// Forward declarations for cross-module SWIG type resolution
// ============================================================================
// These enable typemaps for OpenEye types whose definitions live in the
// OpenEye SWIG runtime (v4). Only types you actually use in your wrapped API
// need full #include — forward declarations suffice for the typemaps.

namespace OEChem {
    class OEMolBase;
    class OEMCMolBase;
    class OEMol;
    class OEGraphMol;
    class OEAtomBase;
    class OEBondBase;
    class OEConfBase;
    class OEMatchBase;
    class OEMolDatabase;
    class oemolistream;
    class oemolostream;
    class OEQMol;
    class OEResidue;
    class OEUniMolecularRxn;
}

namespace OEBio {
    class OEDesignUnit;
    class OEHierView;
    class OEHierResidue;
    class OEHierFragment;
    class OEHierChain;
    class OEInteractionHint;
    class OEInteractionHintContainer;
}

namespace OEDocking {
    class OEReceptor;
}

namespace OEPlatform {
    class oeifstream;
    class oeofstream;
    class oeisstream;
    class oeosstream;
}

namespace OESystem {
    class OERecord;
    class OEMolRecord;
}

// ============================================================================
// Cross-runtime SWIG compatibility layer
// ============================================================================
// OpenEye's Python bindings use SWIG runtime v4; our module uses v5.
// Since the runtimes are separate, SWIG_TypeQuery cannot access OpenEye types.
// We use Python isinstance for type safety and extract the underlying C++
// pointer through SwigPyObject's public __int__ (PyNumber_Long), which both
// runtimes expose. This avoids depending on the private SwigPyObject struct
// layout: int(obj.this) returns the wrapped void* as a Python integer.
//
// Version contract: validated against OpenEye Toolkits 2025.2 (SWIG v4
// wrappers) consumed by this module's SWIG v5 runtime. If a future OpenEye
// release stops exposing SwigPyObject.__int__, extraction fails closed with a
// clean Python exception rather than dereferencing an unknown layout.
//
// This approach enables passing OpenEye objects between Python and C++ without
// serialization. The macros below generate the boilerplate for each type.

%{
// Extract the wrapped C++ pointer from an OpenEye SWIG object via the public
// SwigPyObject.__int__ conversion. Returns NULL (with no Python error pending)
// on any failure, so the caller's typemap fails closed with its own clean
// type error instead of using a bad pointer.
static void* _oemmpa_extract_swig_ptr(PyObject* obj) {
    PyObject* thisAttr = PyObject_GetAttrString(obj, "this");
    if (!thisAttr) {
        // No `this` attribute: not a SWIG-wrapped object. Clear the lookup
        // error so the caller's typemap raises its own clean type error.
        PyErr_Clear();
        return NULL;
    }
    // int(obj.this) yields the underlying void* for a SwigPyObject. This is
    // SWIG's documented, layout-independent extraction path.
    PyObject* asInt = PyNumber_Long(thisAttr);
    Py_DECREF(thisAttr);
    if (!asInt) {
        PyErr_Clear();
        return NULL;
    }
    void* ptr = PyLong_AsVoidPtr(asInt);
    Py_DECREF(asInt);
    if (ptr == NULL && PyErr_Occurred()) {
        PyErr_Clear();
        return NULL;
    }
    return ptr;
}

// ---- Type checker generator macro ----
// Generates a cached isinstance checker for an OpenEye Python type.
// TAG:    identifier suffix (e.g., oemolbase)
// MODULE: Python module string (e.g., "openeye.oechem")
// CLASS:  Python class name string (e.g., "OEMolBase")
#define DEFINE_OE_TYPE_CHECKER(TAG, MODULE, CLASS) \
    static PyObject* _oemmpa_oe_##TAG##_type = NULL; \
    static bool _oemmpa_is_##TAG(PyObject* obj) { \
        if (!_oemmpa_oe_##TAG##_type) { \
            PyObject* mod = PyImport_ImportModule(MODULE); \
            if (mod) { \
                _oemmpa_oe_##TAG##_type = PyObject_GetAttrString(mod, CLASS); \
                Py_DECREF(mod); \
            } \
            if (!_oemmpa_oe_##TAG##_type) return false; \
        } \
        return PyObject_IsInstance(obj, _oemmpa_oe_##TAG##_type) == 1; \
    }

// ---- Molecule types (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oemolbase,    "openeye.oechem", "OEMolBase")
DEFINE_OE_TYPE_CHECKER(oemcmolbase,  "openeye.oechem", "OEMCMolBase")
DEFINE_OE_TYPE_CHECKER(oemol,        "openeye.oechem", "OEMol")
DEFINE_OE_TYPE_CHECKER(oegraphmol,   "openeye.oechem", "OEGraphMol")
DEFINE_OE_TYPE_CHECKER(oeqmol,       "openeye.oechem", "OEQMol")

// ---- Atom / bond / conformer / residue (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oeatombase,   "openeye.oechem", "OEAtomBase")
DEFINE_OE_TYPE_CHECKER(oebondbase,   "openeye.oechem", "OEBondBase")
DEFINE_OE_TYPE_CHECKER(oeconfbase,   "openeye.oechem", "OEConfBase")
DEFINE_OE_TYPE_CHECKER(oeresidue,    "openeye.oechem", "OEResidue")
DEFINE_OE_TYPE_CHECKER(oematchbase,  "openeye.oechem", "OEMatchBase")

// ---- Molecule I/O (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oemolistream, "openeye.oechem", "oemolistream")
DEFINE_OE_TYPE_CHECKER(oemolostream, "openeye.oechem", "oemolostream")
DEFINE_OE_TYPE_CHECKER(oemoldatabase,"openeye.oechem", "OEMolDatabase")

// ---- Reactions (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oeunimolecularrxn, "openeye.oechem", "OEUniMolecularRxn")

// ---- Platform streams (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oeifstream,   "openeye.oechem", "oeifstream")
DEFINE_OE_TYPE_CHECKER(oeofstream,   "openeye.oechem", "oeofstream")
DEFINE_OE_TYPE_CHECKER(oeisstream,   "openeye.oechem", "oeisstream")
DEFINE_OE_TYPE_CHECKER(oeosstream,   "openeye.oechem", "oeosstream")

// ---- Records (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oerecord,     "openeye.oechem", "OERecord")
DEFINE_OE_TYPE_CHECKER(oemolrecord,  "openeye.oechem", "OEMolRecord")

// ---- Bio / hierarchy (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oedesignunit, "openeye.oechem", "OEDesignUnit")
DEFINE_OE_TYPE_CHECKER(oehierview,   "openeye.oechem", "OEHierView")
DEFINE_OE_TYPE_CHECKER(oehierresidue,"openeye.oechem", "OEHierResidue")
DEFINE_OE_TYPE_CHECKER(oehierfragment,"openeye.oechem","OEHierFragment")
DEFINE_OE_TYPE_CHECKER(oehierchain,  "openeye.oechem", "OEHierChain")
DEFINE_OE_TYPE_CHECKER(oeinteractionhint,          "openeye.oechem", "OEInteractionHint")
DEFINE_OE_TYPE_CHECKER(oeinteractionhintcontainer, "openeye.oechem", "OEInteractionHintContainer")

// ---- Docking (openeye.oedocking) ----
DEFINE_OE_TYPE_CHECKER(oereceptor,   "openeye.oedocking", "OEReceptor")

#undef DEFINE_OE_TYPE_CHECKER
%}

// ============================================================================
// Typemap generator macros
// ============================================================================

// Generate const-ref and non-const-ref typemaps for a cross-runtime OpenEye type.
// CPP_TYPE: fully qualified C++ type (e.g., OEChem::OEMolBase)
// CHECKER:  isinstance checker function name
// ERR_MSG:  error message on type mismatch
%define OE_CROSS_RUNTIME_REF_TYPEMAPS(CPP_TYPE, CHECKER, ERR_MSG)

%typemap(in) const CPP_TYPE& (void *argp = 0, int res = 0) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
    if (!SWIG_IsOK(res)) {
        if (CHECKER($input)) {
            argp = _oemmpa_extract_swig_ptr($input);
            if (argp) res = SWIG_OK;
        }
    }
    if (!SWIG_IsOK(res)) {
        SWIG_exception_fail(SWIG_ArgError(res), ERR_MSG);
    }
    if (!argp) {
        SWIG_exception_fail(SWIG_NullReferenceError, "Null reference.");
    }
    $1 = reinterpret_cast< $1_ltype >(argp);
}

%typemap(typecheck, precedence=10) const CPP_TYPE& {
    void *vptr = 0;
    int res = SWIG_ConvertPtr($input, &vptr, $descriptor, SWIG_POINTER_NO_NULL);
    $1 = SWIG_IsOK(res) ? 1 : CHECKER($input) ? 1 : 0;
}

%typemap(in) CPP_TYPE& (void *argp = 0, int res = 0) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
    if (!SWIG_IsOK(res)) {
        if (CHECKER($input)) {
            argp = _oemmpa_extract_swig_ptr($input);
            if (argp) res = SWIG_OK;
        }
    }
    if (!SWIG_IsOK(res)) {
        SWIG_exception_fail(SWIG_ArgError(res), ERR_MSG);
    }
    if (!argp) {
        SWIG_exception_fail(SWIG_NullReferenceError, "Null reference.");
    }
    $1 = reinterpret_cast< $1_ltype >(argp);
}

%typemap(typecheck, precedence=10) CPP_TYPE& {
    void *vptr = 0;
    int res = SWIG_ConvertPtr($input, &vptr, $descriptor, SWIG_POINTER_NO_NULL);
    $1 = SWIG_IsOK(res) ? 1 : CHECKER($input) ? 1 : 0;
}

%enddef

// Generate nullable-pointer typemaps (accepts None) for a cross-runtime type.
%define OE_CROSS_RUNTIME_NULLABLE_PTR_TYPEMAPS(CPP_TYPE, CHECKER, ERR_MSG)

%typemap(in) const CPP_TYPE* (void *argp = 0, int res = 0) {
    if ($input == Py_None) {
        $1 = NULL;
    } else {
        res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
        if (!SWIG_IsOK(res)) {
            if (CHECKER($input)) {
                argp = _oemmpa_extract_swig_ptr($input);
                if (argp) res = SWIG_OK;
            }
        }
        if (!SWIG_IsOK(res)) {
            SWIG_exception_fail(SWIG_ArgError(res), ERR_MSG);
        }
        $1 = reinterpret_cast< $1_ltype >(argp);
    }
}

%typemap(typecheck, precedence=10) const CPP_TYPE* {
    if ($input == Py_None) {
        $1 = 1;
    } else {
        void *vptr = 0;
        int res = SWIG_ConvertPtr($input, &vptr, $descriptor, 0);
        $1 = SWIG_IsOK(res) ? 1 : CHECKER($input) ? 1 : 0;
    }
}

%enddef

// ============================================================================
// Typemap declarations for all OpenEye types
// ============================================================================
// Each type gets const-ref and non-const-ref typemaps. Types that commonly
// appear as optional parameters also get nullable-pointer typemaps.
// These are inert until a wrapped function signature uses the type.

// ---- Molecule hierarchy (OEChem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMolBase,    _oemmpa_is_oemolbase,    "Expected OEMolBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMCMolBase,  _oemmpa_is_oemcmolbase,  "Expected OEMCMolBase-derived object (OEMCMolBase or OEMol).")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMol,        _oemmpa_is_oemol,        "Expected OEMol object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEGraphMol,   _oemmpa_is_oegraphmol,   "Expected OEGraphMol object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEQMol,       _oemmpa_is_oeqmol,       "Expected OEQMol object.")

// ---- Atom / bond / conformer / residue / match (OEChem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEAtomBase,   _oemmpa_is_oeatombase,   "Expected OEAtomBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEBondBase,   _oemmpa_is_oebondbase,   "Expected OEBondBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEConfBase,   _oemmpa_is_oeconfbase,   "Expected OEConfBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEResidue,    _oemmpa_is_oeresidue,    "Expected OEResidue object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMatchBase,  _oemmpa_is_oematchbase,  "Expected OEMatchBase-derived object.")

// ---- Molecule I/O (OEChem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::oemolistream,  _oemmpa_is_oemolistream, "Expected oemolistream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::oemolostream,  _oemmpa_is_oemolostream, "Expected oemolostream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMolDatabase, _oemmpa_is_oemoldatabase,"Expected OEMolDatabase object.")

// ---- Reactions (OEChem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEUniMolecularRxn, _oemmpa_is_oeunimolecularrxn, "Expected OEUniMolecularRxn object.")

// ---- Platform streams (OEPlatform) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeifstream, _oemmpa_is_oeifstream, "Expected oeifstream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeofstream, _oemmpa_is_oeofstream, "Expected oeofstream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeisstream, _oemmpa_is_oeisstream, "Expected oeisstream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeosstream, _oemmpa_is_oeosstream, "Expected oeosstream object.")

// ---- Records (OESystem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OESystem::OERecord,    _oemmpa_is_oerecord,    "Expected OERecord object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OESystem::OEMolRecord, _oemmpa_is_oemolrecord, "Expected OEMolRecord object.")

// ---- Bio / hierarchy (OEBio) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEDesignUnit,   _oemmpa_is_oedesignunit, "Expected OEDesignUnit object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierView,     _oemmpa_is_oehierview,   "Expected OEHierView object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierResidue,  _oemmpa_is_oehierresidue,"Expected OEHierResidue object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierFragment,  _oemmpa_is_oehierfragment,"Expected OEHierFragment object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierChain,    _oemmpa_is_oehierchain,  "Expected OEHierChain object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEInteractionHint,          _oemmpa_is_oeinteractionhint,          "Expected OEInteractionHint object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEInteractionHintContainer, _oemmpa_is_oeinteractionhintcontainer, "Expected OEInteractionHintContainer object.")

// ---- Docking (OEDocking) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEDocking::OEReceptor, _oemmpa_is_oereceptor, "Expected OEReceptor object.")

// ============================================================================
// Include STL typemaps
// ============================================================================
%include "std_string.i"
%include "std_vector.i"
%include "stdint.i"
%include "exception.i"

// ============================================================================
// Exception handling
// ============================================================================
%exception {
    try {
        $action
    } catch (const OEMMPA::OEMMPAError& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (...) {
        SWIG_exception(SWIG_RuntimeError, "Unknown C++ exception");
    }
}

// ============================================================================
// Template instantiations for container types
// ============================================================================
%template(StringVector) std::vector<std::string>;
%template(UnsignedIntVector) std::vector<unsigned int>;
%template(UnsignedIntVectorVector) std::vector<std::vector<unsigned int>>;
%template(LoadErrorVector) std::vector<OEMMPA::LoadError>;
%template(AcceptedMoleculeVector) std::vector<OEMMPA::AcceptedMolecule>;
%template(CutBondVector) std::vector<OEMMPA::CutBond>;
%template(FragmentationVector) std::vector<OEMMPA::Fragmentation>;
%template(EnvironmentFingerprintVector) std::vector<OEMMPA::EnvironmentFingerprint>;
%template(QueryEnvironmentVector) std::vector<OEMMPA::QueryEnvironment>;
%template(PairEnvironmentSmirksVector) std::vector<OEMMPA::PairEnvironmentSmirks>;
%template(MatchedPairVector) std::vector<OEMMPA::MatchedPair>;
%template(RuleEnvironmentStatisticsVector) std::vector<OEMMPA::RuleEnvironmentStatistics>;
%template(TransformVector) std::vector<OEMMPA::Transform>;
%template(TransformProductVector) std::vector<OEMMPA::TransformProduct>;
%template(GeneratedProductVector) std::vector<OEMMPA::GeneratedProduct>;
// A series and a member carry their contents from construction, so neither has a
// default constructor. This is the SWIG library's opt-out for such an element
// type: it drops the two vector members that would have to default-construct
// one, the sized constructor vector(size_type) and the one-argument
// resize(size_type). The two-argument resize(size_type, const value_type&) is
// unaffected and stays wrapped.
%std_nodefconst_type(OEMMPA::SeriesMember);
%std_nodefconst_type(OEMMPA::Series);
%template(SeriesMemberVector) std::vector<OEMMPA::SeriesMember>;
%template(SeriesVector) std::vector<OEMMPA::Series>;
// The query records are plain aggregates with default constructors, so the
// stock vector wrappers apply; std::vector<double> is new to this module.
%template(DoubleVector) std::vector<double>;
%template(SubstituentRankingVector) std::vector<OEMMPA::SubstituentRanking>;
%template(SubstituentRankingEvidenceVector) std::vector<OEMMPA::SubstituentRankingEvidence>;
%template(SeriesSimilarityVector) std::vector<OEMMPA::SeriesSimilarity>;
%template(SeriesSimilarityEvidenceVector) std::vector<OEMMPA::SeriesSimilarityEvidence>;

// ============================================================================
// Version macros
// ============================================================================
#define OEMMPA_VERSION_MAJOR 5
#define OEMMPA_VERSION_MINOR 0
#define OEMMPA_VERSION_PATCH 0

// ============================================================================
// GIL release opt-ins (MUST precede the %include of Analyzer.h / DuckDBStore.h)
// ============================================================================
// These 10 methods are the ONLY wrappers that release the GIL.
%threadallow OEMMPA::Analyzer::Analyze;
%threadallow OEMMPA::Analyzer::GetPairs;
%threadallow OEMMPA::Analyzer::GetSeries;
%threadallow OEMMPA::Analyzer::GetTransforms;
%threadallow OEMMPA::Analyzer::SaveTo;
%threadallow OEMMPA::DuckDBStore::AddMoleculesFromSmilesFile;
%threadallow OEMMPA::DuckDBStore::AddPropertiesFromCsvFile;
%threadallow OEMMPA::DuckDBStore::GetSeries;
// Both take their request by value so the wrapper copies it before releasing
// the GIL (the plan 1 ruling for GetSeries(SeriesSelection)).
%threadallow OEMMPA::DuckDBStore::RankNextSubstituents;
%threadallow OEMMPA::DuckDBStore::FindSimilarSeries;

// ============================================================================
// Phase 1 public headers
// ============================================================================
%ignore OEMMPA::FragmentationStrategy::Clone;
%ignore OEMMPA::SmartsFragmentationStrategy::Clone;
%ignore OEMMPA::BondIndexFragmentationStrategy::Clone;
%ignore OEMMPA::MoleculeRecord::GetMol;
%ignore OEMMPA::AnalysisMethod::GetFragmenter;
%ignore OEMMPA::AnalysisMethod::SetFragmenter;
%ignore OEMMPA::FragmentationMethod::GetFragmenter;
%ignore OEMMPA::FragmentationMethod::SetFragmenter;
// build_transforms is an implementation detail shared by the C++ query paths;
// exporting it would add a public Python name this binding never intended.
%ignore OEMMPA::build_transforms;
// Python has no rvalue concept, so SWIG would wrap the rvalue overload as a
// second AddPair indistinguishable from the const-reference one and dispatch
// between them. Ignore it outright so Python keeps reaching the copying
// overload.
%ignore OEMMPA::Transform::AddPair(MatchedPair &&);

%include "../include/oemmpa/Error.h"
%include "../include/oemmpa/DatabaseSummary.h"
%include "../include/oemmpa/EnvironmentFingerprint.h"
// Desalting lives in the standalone oedesalt library. oemmpa wraps its
// OEDESALT::Desalter into this single SWIG module so the facade can construct a
// desalter and thread it through the C++ operations. SaltPattern/DesaltResult
// hold OpenEye types with no cross-runtime typemap, and Desalt() returns an
// OEGraphMol by value that cannot be handed back safely, so only the factories
// (WithBundledPatterns / FromFiles), PatternCount, and IsAggressive are exposed.
// The raw vector<SaltPattern> constructor is unusable from Python and ignored.
%ignore OEDESALT::SaltPattern;
%ignore OEDESALT::DesaltResult;
%ignore OEDESALT::load_salt_patterns;
%ignore OEDESALT::Desalter::Desalt;
%ignore OEDESALT::Desalter::Desalter;
%include "oedesalt/Desalter.h"
%include "../include/oemmpa/QueryEnvironment.h"
%include "../include/oemmpa/Fragmentation.h"
%include "../include/oemmpa/MoleculeRecord.h"
%include "../include/oemmpa/LoadReport.h"
// kHydrogenVariableSmiles is the shared C++ literal for the synthesized hydrogen
// variable fragment. SWIG turns a namespace-scope global into a module attribute
// on the generated proxy plus an entry on the `cvar` holder; no Python caller
// wants either. The directive has to sit here rather than with the FragmentSmiles
// ignores below, because a %ignore placed after the %include of the header that
// declares the symbol has already missed its chance.
%ignore OEMMPA::kHydrogenVariableSmiles;
// count_variable_heavy_atoms is the C++ seam that keeps the store, the
// in-memory index and MatchedPair's accessors counting a fragment the same way.
// Its std::optional return has no typemap, so SWIG would expose an opaque proxy
// object; Python callers read the counts off the pair instead. Same placement
// rule as the ignore above: it must precede the %include.
%ignore OEMMPA::count_variable_heavy_atoms;
%include "../include/oemmpa/MatchedPair.h"
%include "../include/oemmpa/Transform.h"
%include "../include/oemmpa/RuleEnvironmentStatistics.h"
%include "../include/oemmpa/TransformApplication.h"
%include "../include/oemmpa/QueryOptions.h"
// Only the fragment renderer, so a Python measurement can compare OEMMPA's
// stereo cleanup against RDKit's parse-time cleanup on the same molecule.
// Everything else in the header is C++-only plumbing: StereoCleanupContext and
// TetrahedralParity hold raw OEAtomBase pointers and std::vector<bool> with no
// cross-runtime typemap, assert_tetrahedral_parity() takes a std::function, the
// three predicates are meaningless without a context to pass them, and
// ScopedMapIndexClear is a C++ scope guard whose whole contract is the point at
// which its destructor runs. clear_fragmentation_invalidated_stereo() is the
// cleanup pass render_fragment_smiles() already runs for its caller; exporting
// it as well would offer Python a second, in-place entry point into the header
// and move the fragment SMILES flag set out of the one function that owns it.
%ignore OEMMPA::StereoCleanupContext;
%ignore OEMMPA::TetrahedralParity;
%ignore OEMMPA::ScopedMapIndexClear;
%ignore OEMMPA::build_stereo_cleanup_context;
%ignore OEMMPA::rdkit_keeps_tetrahedral_tag;
%ignore OEMMPA::has_dependent_ring_stereo;
%ignore OEMMPA::has_specified_stereo;
%ignore OEMMPA::is_potential_stereocenter;
%ignore OEMMPA::capture_tetrahedral_parity;
%ignore OEMMPA::assert_tetrahedral_parity;
%ignore OEMMPA::clear_fragmentation_invalidated_stereo;
// Internal ordering helper shared by the fragmenter and the query canonicalisation; not part of the binding.
%ignore OEMMPA::order_atoms_for_label_stable_rendering;
%include "../include/oemmpa/FragmentSmiles.h"
// The ported ring perception, so a Python measurement can score OEMMPA's
// symmetrized SSSR and bridgehead set against RDKit's Chem.GetSymmSSSR() and
// rdqueries.IsBridgeheadQueryAtom() on the same molecule. find_bridgeheads()
// returns a std::vector<bool> with no cross-runtime typemap and is redundant
// with bridgehead_atoms(), which is the form the measurement consumes.
%ignore OEMMPA::find_bridgeheads;
// The queue-bound seam is a test fixture, not API: it exists so a C++ test can
// reach the abort path without carrying a molecule big enough to queue 200000
// partial rings, and it is process-wide, so a Python caller could only use it to
// break every other thread's ring perception. Same placement rule as the ignores
// above -- it must precede the %include.
%ignore OEMMPA::detail::set_ring_search_queue_bound_for_tests;
%ignore OEMMPA::detail::ring_search_queue_bound;
%include "../include/oemmpa/RingPerception.h"
%include "../include/oemmpa/FragmentationStrategy.h"
%include "../include/oemmpa/Fragmenter.h"
// A member's enumeration labels are an enum class per compound, and this module
// has no vector wrapper for that type: the accessor and the constructor that
// takes the labels would surface in Python as opaque pointers. Same placement
// rule as the ignores above: before the %include.
%ignore OEMMPA::SeriesMember::GetEnumerationLabels;
%ignore OEMMPA::SeriesMember::SeriesMember(std::string, unsigned int, std::vector<unsigned int>, std::vector<std::string>, std::vector<Fragmentation::EnumerationLabel>);
%include "../include/oemmpa/Series.h"
// The store validates its own support ladders, so Python never calls this one.
%ignore OEMMPA::validate_support_ladder;
// std::string members are natural by default, std::vector members are not: SWIG
// wraps them as a view onto the live member, which is wrong in both directions.
// Writing one would only accept another wrapped vector, and reading one hands
// Python a view whose lifetime is the record's: the store returns its results by
// value, so `store.RankNextSubstituents(request).candidates` would outlive the
// record it points into and read freed memory as an empty sequence, which no
// caller can tell apart from "no candidates". Every vector member of these four
// records is therefore natural: a write accepts any Python sequence, and a read
// returns an independent copy that stays valid after the record is gone. The cost
// is one copy per access, which the Python helpers pay once when they materialise
// a result.
%naturalvar OEMMPA::SubstituentRankingRequest::query_variable_smiles;
%naturalvar OEMMPA::SubstituentRankingRequest::query_values;
%naturalvar OEMMPA::SubstituentRankingRequest::support_ladder;
%naturalvar OEMMPA::SubstituentRankingResult::query_used;
%naturalvar OEMMPA::SubstituentRankingResult::dropped;
%naturalvar OEMMPA::SubstituentRankingResult::candidates;
%naturalvar OEMMPA::SubstituentRankingResult::evidence;
%naturalvar OEMMPA::SeriesSimilarityRequest::query_variable_smiles;
%naturalvar OEMMPA::SeriesSimilarityRequest::query_values;
%naturalvar OEMMPA::SeriesSimilarityResult::matches;
%naturalvar OEMMPA::SeriesSimilarityResult::evidence;
// `selection` stays a view onto the live member, so `request.selection.cut_count = 2`
// still reaches the request, but SWIG's getter hands out a bare pointer into the parent
// without retaining it: a view read off a temporary request, as in
// `SubstituentRankingRequest().selection`, outlives the record it points into and then
// reads and writes whatever the allocator gives the next request. These overrides keep
// the generated getter and setter and only pin the parent to the view they hand back,
// which is enough because the view is the only thing that can reach the freed record.
%extend OEMMPA::SubstituentRankingRequest {
%pythoncode %{
    def _selection_view(self):
        r"""Return the series bounds of this request as a live view.

        :returns: The request's own :class:`SeriesSelection`; writes through it change
            this request, which the view keeps alive for as long as it exists.
        """
        view = _oemmpa.SubstituentRankingRequest_selection_get(self)
        view._owner = self
        return view

    selection = property(_selection_view, _oemmpa.SubstituentRankingRequest_selection_set,
                         doc="Series bounds, a live view into this request (spec 7.5).")
%}
}
%extend OEMMPA::SeriesSimilarityRequest {
%pythoncode %{
    def _selection_view(self):
        r"""Return the series bounds of this request as a live view.

        :returns: The request's own :class:`SeriesSelection`; writes through it change
            this request, which the view keeps alive for as long as it exists.
        """
        view = _oemmpa.SeriesSimilarityRequest_selection_get(self)
        view._owner = self
        return view

    selection = property(_selection_view, _oemmpa.SeriesSimilarityRequest_selection_set,
                         doc="Series bounds, a live view into this request (spec 7.5).")
%}
}
%include "../include/oemmpa/SeriesQuery.h"
// Internal holder for the index's cache lock; never part of the binding surface (it would hand out a raw std::mutex).
%ignore OEMMPA::CacheMutex;
%include "../include/oemmpa/MemoryIndex.h"
%include "../include/oemmpa/PairScoring.h"
%include "../include/oemmpa/Statistics.h"
%include "../include/oemmpa/AnalysisMethod.h"
%include "../include/oemmpa/FragmentationMethod.h"
#if OEMMPA_HAS_DUCKDB
%include "../include/oemmpa/DuckDBStore.h"
#endif
%include "../include/oemmpa/Analyzer.h"
%include "../include/oemmpa/oemmpa.h"

// ============================================================================
// Bulk result export (single-crossing to_dicts)
// ============================================================================
// Materializing pairs/transforms one field at a time costs ~10 Python<->C++
// getter crossings per row, which dominates to_dicts() at scale. These helpers
// build the entire list of dicts in one C call. The implementations live in a
// plain %{ %} block (compiled into the wrapper, not itself wrapped); only the
// two public functions below are exposed to Python.
%{
// value is a new reference (or NULL on allocation failure); on success the dict
// takes its own reference, so we always release ours. Returns -1 on any failure
// with our reference dropped, so callers can bail without leaking.
static int _oemmpa_dict_set_new(PyObject* dict, const char* key, PyObject* value) {
    if (value == NULL) {
        return -1;
    }
    const int rc = PyDict_SetItemString(dict, key, value);
    Py_DECREF(value);
    return rc;
}

// Mirror PairResult._id_or_internal: the external id string when non-empty,
// otherwise the internal integer id.
static PyObject* _oemmpa_pair_id_object(const std::string& external_id, unsigned int internal_id) {
    if (external_id.empty()) {
        return PyLong_FromUnsignedLong(internal_id);
    }
    return PyUnicode_FromStringAndSize(
        external_id.data(), static_cast<Py_ssize_t>(external_id.size()));
}

static PyObject* _oemmpa_str_object(const std::string& value) {
    return PyUnicode_FromStringAndSize(
        value.data(), static_cast<Py_ssize_t>(value.size()));
}

// Byte-for-byte equivalent to [PairResult(p).to_dict() for p in pairs].
static PyObject* matched_pair_dicts(const std::vector<OEMMPA::MatchedPair>& pairs) {
    PyObject* list = PyList_New(static_cast<Py_ssize_t>(pairs.size()));
    if (list == NULL) {
        return NULL;
    }
    Py_ssize_t index = 0;
    for (const OEMMPA::MatchedPair& pair : pairs) {
        PyObject* row = PyDict_New();
        if (row == NULL) {
            Py_DECREF(list);
            return NULL;
        }

        // Populate all scalar fields first (so environment_smirks cleanup is isolated)
        if (_oemmpa_dict_set_new(row, "source_id",
                _oemmpa_pair_id_object(pair.GetSourceExternalId(), pair.GetSourceMoleculeId())) < 0
            || _oemmpa_dict_set_new(row, "target_id",
                _oemmpa_pair_id_object(pair.GetTargetExternalId(), pair.GetTargetMoleculeId())) < 0
            || _oemmpa_dict_set_new(row, "constant",
                _oemmpa_str_object(pair.GetConstantSmiles())) < 0
            || _oemmpa_dict_set_new(row, "source_variable",
                _oemmpa_str_object(pair.GetSourceVariableSmiles())) < 0
            || _oemmpa_dict_set_new(row, "target_variable",
                _oemmpa_str_object(pair.GetTargetVariableSmiles())) < 0
            || _oemmpa_dict_set_new(row, "transform",
                _oemmpa_str_object(pair.GetTransformSmiles())) < 0
            || _oemmpa_dict_set_new(row, "cut_count",
                PyLong_FromUnsignedLong(pair.GetCutCount())) < 0
            || _oemmpa_dict_set_new(row, "heavy_atom_delta",
                PyLong_FromLong(pair.GetHeavyAtomDelta())) < 0
            || _oemmpa_dict_set_new(row, "heavy_bond_delta",
                PyLong_FromLong(pair.GetHeavyBondDelta())) < 0) {
            Py_DECREF(row);
            Py_DECREF(list);
            return NULL;
        }

        // Build environment_smirks list (immediately before transfer to avoid intermediate error paths)
        const auto& env_smirks = pair.GetEnvironmentSmirks();
        PyObject* smirks_list = PyList_New(static_cast<Py_ssize_t>(env_smirks.size()));
        if (smirks_list == NULL) {
            Py_DECREF(row);
            Py_DECREF(list);
            return NULL;
        }
        for (size_t i = 0; i < env_smirks.size(); ++i) {
            PyObject* entry = PyDict_New();
            if (entry == NULL
                || _oemmpa_dict_set_new(entry, "radius",
                    PyLong_FromUnsignedLong(env_smirks[i].radius)) < 0
                || _oemmpa_dict_set_new(entry, "smirks",
                    _oemmpa_str_object(env_smirks[i].smirks)) < 0) {
                Py_XDECREF(entry);
                Py_DECREF(smirks_list);
                Py_DECREF(row);
                Py_DECREF(list);
                return NULL;
            }
            PyList_SetItem(smirks_list, i, entry);  // steals the entry reference
        }

        // Transfer smirks_list to row (steals reference; self-cleans on failure; no further decref needed)
        if (_oemmpa_dict_set_new(row, "environment_smirks", smirks_list) < 0) {
            Py_DECREF(row);
            Py_DECREF(list);
            return NULL;
        }
        PyList_SetItem(list, index, row);  // steals the row reference
        ++index;
    }
    return list;
}

// Byte-for-byte equivalent to the Python series_rows(): one dict per
// (member, compound) in the order of the vector, of each series' members and of
// their compound ids. Nine getter crossings per row is what this replaces, and
// reading the C++ vectors directly also avoids the non-owning SWIG views that
// GetMembers() and GetMoleculeIds() hand back. The vector arrives by copy from a
// Python sequence, so callers pass only the series they need converted (the
// notebook preview passes one at a time). Failure anywhere releases everything
// built so far and returns NULL with the CPython error set, never a half-built
// list: Py_DECREF of a list whose tail slots are still NULL is safe, because
// list deallocation skips NULL items, which is why the slots are filled in order.
static PyObject* series_row_dicts(const std::vector<OEMMPA::Series>& series) {
    Py_ssize_t total = 0;
    for (const OEMMPA::Series& entry : series) {
        for (const OEMMPA::SeriesMember& member : entry.GetMembers()) {
            total += static_cast<Py_ssize_t>(member.GetMoleculeIds().size());
        }
    }
    PyObject* list = PyList_New(total);
    if (list == NULL) {
        return NULL;
    }
    // A member whose external ids are shorter than its molecule ids is
    // malformed; this keeps the read in bounds and yields the empty string the
    // Python row builder also produced. Held by const reference so neither arm
    // of the ternary below copies.
    const std::string missing_external_id;
    Py_ssize_t index = 0;
    for (const OEMMPA::Series& entry : series) {
        for (const OEMMPA::SeriesMember& member : entry.GetMembers()) {
            const std::vector<unsigned int>& ids = member.GetMoleculeIds();
            const std::vector<std::string>& external_ids = member.GetExternalIds();
            for (size_t i = 0; i < ids.size(); ++i) {
                PyObject* row = PyDict_New();
                if (row == NULL) {
                    Py_DECREF(list);
                    return NULL;
                }
                // Every value is built inside the chain so that a short-circuit
                // on an earlier failure cannot strand an unowned new reference.
                if (_oemmpa_dict_set_new(row, "series_id",
                        PyLong_FromUnsignedLongLong(entry.GetSeriesId())) < 0
                    || _oemmpa_dict_set_new(row, "constant",
                        _oemmpa_str_object(entry.GetConstantSmiles())) < 0
                    || _oemmpa_dict_set_new(row, "cut_count",
                        PyLong_FromUnsignedLong(entry.GetCutCount())) < 0
                    || _oemmpa_dict_set_new(row, "constant_num_heavies",
                        PyLong_FromUnsignedLong(entry.GetConstantHeavyAtomCount())) < 0
                    || _oemmpa_dict_set_new(row, "num_members",
                        PyLong_FromSize_t(entry.Size())) < 0
                    || _oemmpa_dict_set_new(row, "rgroup",
                        _oemmpa_str_object(member.GetVariableSmiles())) < 0
                    || _oemmpa_dict_set_new(row, "rgroup_num_heavies",
                        PyLong_FromUnsignedLong(member.GetHeavyAtomCount())) < 0
                    || _oemmpa_dict_set_new(row, "molecule_id",
                        PyLong_FromUnsignedLong(ids[i])) < 0
                    || _oemmpa_dict_set_new(row, "compound_id",
                        _oemmpa_str_object(i < external_ids.size()
                                           ? external_ids[i]
                                           : missing_external_id)) < 0) {
                    Py_DECREF(row);
                    Py_DECREF(list);
                    return NULL;
                }
                PyList_SetItem(list, index, row);  // steals the row reference
                ++index;
            }
        }
    }
    return list;
}
%}

// Expose only the public builders (SWIG passes the PyObject* return through).
// Transforms are intentionally not bulk-exported: with only two fields the
// per-row getter count is too low to beat the C dict-building overhead.
PyObject* matched_pair_dicts(const std::vector<OEMMPA::MatchedPair>& pairs);
PyObject* series_row_dicts(const std::vector<OEMMPA::Series>& series);

// ============================================================================
// Module-level Python convenience code
// ============================================================================
%pythoncode %{
# Derive the proxy module's __version__ from the compiled OEMMPA_VERSION_*
# macros (wrapped from oemmpa.h) so it never drifts from the released version.
__version__ = "%d.%d.%d" % (OEMMPA_VERSION_MAJOR, OEMMPA_VERSION_MINOR, OEMMPA_VERSION_PATCH)
%}
