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
    # wrapper bodies; every absence assertion uses a tighter window so it cannot
    # bleed into the following wrapper and read a neighbor's macro. The two query
    # wrappers convert a request of ten fields before they release, which puts
    # their macro ~1500 chars in, still far short of the next wrapper's.
    RELEASE_SCAN = 1800
    GETTER_SCAN = 800
    def wrapper_releases(pattern, scan=RELEASE_SCAN):
        blocks = re.findall(pattern + r"[\s\S]{0," + str(scan) + r"}", text)
        return blocks and any(any(r in b for r in RELEASE) for b in blocks)
    # Every one of the ten target methods must release the GIL. The two series reads
    # are named by their overload wrappers rather than by the bare method name: the
    # prefix "DuckDBStore_GetSeries" also matches the GIL-held GetSeriesCount wrapper,
    # so a bare prefix would pass on a build where only one of the three releases.
    for name in ("Analyzer_Analyze", "Analyzer_GetPairs", "Analyzer_GetSeries",
                 "Analyzer_GetTransforms", "Analyzer_SaveTo",
                 "DuckDBStore_AddMoleculesFromSmilesFile",
                 "DuckDBStore_AddPropertiesFromCsvFile",
                 "DuckDBStore_GetSeries__SWIG_0", "DuckDBStore_GetSeries__SWIG_1",
                 "DuckDBStore_RankNextSubstituents", "DuckDBStore_FindSimilarSeries"):
        assert wrapper_releases(r"_wrap_" + name), f"{name} should release the GIL"
    # A sampled trivial getter must NOT release (proves selectivity, no churn).
    getter = re.findall(r"_wrap_Analyzer_GetMethodName[\s\S]{0," + str(GETTER_SCAN) + r"}", text)
    assert getter, "expected a getter wrapper to sample"
    assert not any(any(r in b for r in RELEASE) for b in getter), \
        "trivial getter should stay GIL-held"
    # Exclusivity: sample other non-target wrappers and assert they do NOT release.
    # Every name here must exist in the wrapper, otherwise the check is vacuous.
    non_targets = ("Analyzer_GetStrippedNames", "DuckDBStore_GetRowCount",
                   "AnalysisMethod_LastAnalyzeWorkerCount", "Fragmenter_GetMaxCuts",
                   "Analyzer_ProducesSeries", "DuckDBStore_GetSeriesCount")
    for name in non_targets:
        assert re.search(r"_wrap_" + name + r"\b", text), f"{name} is not a wrapper any more"
        assert not wrapper_releases(r"_wrap_" + name, scan=GETTER_SCAN), f"{name} should NOT release the GIL"
    # Total GIL macro count: must be low. Only the BEGIN macros are counted, one
    # per releasing wrapper body, so the bound is the 10 targets times their
    # overloads rather than a begin-and-end pair count.
    # Pin the observed count so a regression that goes broad will fail.
    macro_count = text.count("SWIG_PYTHON_THREAD_BEGIN_ALLOW") + text.count("Py_BEGIN_ALLOW_THREADS")
    assert macro_count <= 24, f"Expected ≤24 GIL macros (10 targets), found {macro_count}"


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


@pytest.mark.skipif(not WRAP.exists(), reason="generated SWIG wrapper not built")
def test_query_requests_are_copied_before_the_gil_is_released():
    # Same hazard as the selection: the request is a Python-owned aggregate whose
    # vectors and bounds are writable, so the GIL-free query must run on a copy.
    import re
    text = WRAP.read_text()
    for wrapper, record in (
        ("DuckDBStore_RankNextSubstituents", "SubstituentRankingRequest"),
        ("DuckDBStore_FindSimilarSeries", "SeriesSimilarityRequest"),
    ):
        source = _wrapper_source(text, wrapper)
        before, macro, after = source.partition("SWIG_PYTHON_THREAD_BEGIN_ALLOW")
        assert macro, f"{wrapper} should release the GIL"
        assert re.search(r"OEMMPA::" + record + r" arg\d+ ;", before), \
            f"{wrapper} must take the request by value"
        assert re.search(r"arg\d+ = \*\w+;", before), \
            f"{wrapper} must copy the request before the GIL is released"
        assert record not in after, \
            f"{wrapper} must not read the Python-owned request after releasing the GIL"


@pytest.mark.skipif(not WRAP.exists(), reason="generated SWIG wrapper not built")
def test_query_options_are_copied_before_the_gil_is_released():
    # Same hazard as the selection and requests: the options are a Python-owned
    # aggregate whose bounds are writable, so the GIL-free call must run on a copy.
    import re
    text = WRAP.read_text()
    for wrapper in ("Analyzer_GetPairs__SWIG_1", "Analyzer_GetTransforms__SWIG_1",
                    "Analyzer_SaveTo__SWIG_1"):
        source = _wrapper_source(text, wrapper)
        before, macro, after = source.partition("SWIG_PYTHON_THREAD_BEGIN_ALLOW")
        assert macro, f"{wrapper} should release the GIL"
        assert re.search(r"OEMMPA::QueryOptions arg\d+ ;", before), \
            f"{wrapper} must take the options by value"
        assert re.search(r"arg\d+ = \*\w+;", before), \
            f"{wrapper} must copy the options before the GIL is released"
        assert "QueryOptions" not in after, \
            f"{wrapper} must not read the Python-owned options after releasing the GIL"
