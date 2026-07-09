import pathlib
import pytest

WRAP = pathlib.Path("build-debug/swig/oemmpaPYTHON_wrap.cxx")

@pytest.mark.skipif(not WRAP.exists(), reason="generated SWIG wrapper not built")
def test_target_methods_release_gil():
    text = WRAP.read_text()
    # The long calls must be bracketed by GIL release.
    assert "SWIG_PYTHON_THREAD_BEGIN_ALLOW" in text or "Py_BEGIN_ALLOW_THREADS" in text
    # Spot-check: the long-running methods release the GIL; a trivial getter does not.
    import re
    RELEASE = ("SWIG_PYTHON_THREAD_BEGIN_ALLOW", "Py_BEGIN_ALLOW_THREADS")
    # Scan windows are the number of chars examined after a wrapper name. The
    # release check uses a wider window so it reaches the GIL macro inside larger
    # wrapper bodies; the getter check (absence assertion) uses a tighter window
    # so it cannot bleed into the following wrapper and read a neighbor's macro.
    RELEASE_SCAN = 1500
    GETTER_SCAN = 800
    def wrapper_releases(pattern):
        blocks = re.findall(pattern + r"[\s\S]{0," + str(RELEASE_SCAN) + r"}", text)
        return blocks and any(any(r in b for r in RELEASE) for b in blocks)
    # Every one of the six target methods must release the GIL.
    for name in ("Analyzer_Analyze", "Analyzer_GetPairs", "Analyzer_GetTransforms",
                 "Analyzer_SaveTo", "DuckDBStore_AddMoleculesFromSmilesFile",
                 "DuckDBStore_AddPropertiesFromCsvFile"):
        assert wrapper_releases(r"_wrap_" + name), f"{name} should release the GIL"
    # A sampled trivial getter must NOT release (proves selectivity, no churn).
    getter = re.findall(r"_wrap_Analyzer_GetMethodName[\s\S]{0," + str(GETTER_SCAN) + r"}", text)
    assert getter, "expected a getter wrapper to sample"
    assert not any(any(r in b for r in RELEASE) for b in getter), \
        "trivial getter should stay GIL-held"
    # Exclusivity: sample other non-target wrappers and assert they do NOT release.
    non_targets = ("Analyzer_GetFragmentationCount", "DuckDBStore_GetDatabasePath",
                   "AnalysisMethod_GetMethodName", "Fragmenter_GetMaxCuts")
    for name in non_targets:
        assert not wrapper_releases(r"_wrap_" + name), f"{name} should NOT release the GIL"
    # Total GIL macro count: must be low (only the 6 targets × overloads, begin+end).
    # Pin the observed count so a regression that goes broad will fail.
    macro_count = text.count("SWIG_PYTHON_THREAD_BEGIN_ALLOW") + text.count("Py_BEGIN_ALLOW_THREADS")
    assert macro_count <= 24, f"Expected ≤24 GIL macros (6 targets), found {macro_count}"
