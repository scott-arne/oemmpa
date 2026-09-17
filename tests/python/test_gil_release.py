import pathlib
import pytest

WRAP = pathlib.Path("build-debug/swig/oemmpaPYTHON_wrap.cxx")


def _wrapper_source(text, name):
    """Return the generated source of a single wrapper function.

    :param text: The whole generated wrapper translation unit.
    :param name: Wrapper function name without its ``_wrap_`` prefix.
    :returns: The text from that wrapper's definition to the next definition.
    """
    start = text.index("SWIGINTERN PyObject *_wrap_" + name + "(")
    return text[start:text.index("\nSWIGINTERN ", start + 1)]


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
    # Every one of the eight target methods must release the GIL. The two series reads
    # are named by their overload wrappers rather than by the bare method name: the
    # prefix "DuckDBStore_GetSeries" also matches the GIL-held GetSeriesCount wrapper,
    # so a bare prefix would pass on a build where only one of the three releases.
    for name in ("Analyzer_Analyze", "Analyzer_GetPairs", "Analyzer_GetSeries",
                 "Analyzer_GetTransforms", "Analyzer_SaveTo",
                 "DuckDBStore_AddMoleculesFromSmilesFile",
                 "DuckDBStore_AddPropertiesFromCsvFile",
                 "DuckDBStore_GetSeries__SWIG_0", "DuckDBStore_GetSeries__SWIG_1"):
        assert wrapper_releases(r"_wrap_" + name), f"{name} should release the GIL"
    # A sampled trivial getter must NOT release (proves selectivity, no churn).
    getter = re.findall(r"_wrap_Analyzer_GetMethodName[\s\S]{0," + str(GETTER_SCAN) + r"}", text)
    assert getter, "expected a getter wrapper to sample"
    assert not any(any(r in b for r in RELEASE) for b in getter), \
        "trivial getter should stay GIL-held"
    # Exclusivity: sample other non-target wrappers and assert they do NOT release.
    non_targets = ("Analyzer_GetFragmentationCount", "DuckDBStore_GetDatabasePath",
                   "AnalysisMethod_GetMethodName", "Fragmenter_GetMaxCuts",
                   "Analyzer_ProducesSeries", "DuckDBStore_GetSeriesCount")
    for name in non_targets:
        assert not wrapper_releases(r"_wrap_" + name), f"{name} should NOT release the GIL"
    # Total GIL macro count: must be low (only the 8 targets × overloads, begin+end).
    # Pin the observed count so a regression that goes broad will fail.
    macro_count = text.count("SWIG_PYTHON_THREAD_BEGIN_ALLOW") + text.count("Py_BEGIN_ALLOW_THREADS")
    assert macro_count <= 24, f"Expected ≤24 GIL macros (8 targets), found {macro_count}"


@pytest.mark.skipif(not WRAP.exists(), reason="generated SWIG wrapper not built")
def test_series_selection_is_copied_before_the_gil_is_released():
    # The selection is a Python object whose four bounds are writable attributes, so
    # the GIL-free read must work from a copy taken while the GIL is still held.
    # Reading it through a pointer after the release would let another thread change
    # the bounds mid-query: DuckDBStore::GetSeries reads three of them when it binds
    # the SQL parameters and min_series_length again for every grouped series, so a
    # torn read could combine bounds that were never one coherent selection.
    import re
    text = WRAP.read_text()

    with_selection = _wrapper_source(text, "DuckDBStore_GetSeries__SWIG_0")
    before, macro, after = with_selection.partition("SWIG_PYTHON_THREAD_BEGIN_ALLOW")
    assert macro, "the selection overload should release the GIL"
    assert re.search(r"OEMMPA::SeriesSelection arg\d+ ;", before), \
        "the selection must reach the call by value, not as a pointer to the Python object"
    assert re.search(r"arg\d+ = \*\w+;", before), \
        "the copy of the selection must be taken before the GIL is released"
    assert "SeriesSelection" not in after, \
        "the GIL-free call must not read the Python-owned selection through a pointer"

    defaulted = _wrapper_source(text, "DuckDBStore_GetSeries__SWIG_1")
    before, macro, after = defaulted.partition("SWIG_PYTHON_THREAD_BEGIN_ALLOW")
    assert macro, "the defaulted overload should release the GIL"
    assert "SeriesSelection" not in after, \
        "the defaulted overload must not touch a Python-owned selection either"
