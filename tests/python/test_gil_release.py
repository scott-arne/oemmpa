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
    def wrapper_releases(pattern):
        blocks = re.findall(pattern + r"[\s\S]{0,1500}", text)
        return blocks and any(any(r in b for r in RELEASE) for b in blocks)
    # Every one of the six target methods must release the GIL.
    for name in ("Analyzer_Analyze", "Analyzer_GetPairs", "Analyzer_GetTransforms",
                 "Analyzer_SaveTo", "DuckDBStore_AddMoleculesFromSmilesFile",
                 "DuckDBStore_AddPropertiesFromCsvFile"):
        assert wrapper_releases(r"_wrap_" + name), f"{name} should release the GIL"
    # A sampled trivial getter must NOT release (proves selectivity, no churn).
    getter = re.findall(r"_wrap_Analyzer_GetMethodName[\s\S]{0,800}", text)
    assert getter, "expected a getter wrapper to sample"
    assert not any(any(r in b for r in RELEASE) for b in getter), \
        "trivial getter should stay GIL-held"
