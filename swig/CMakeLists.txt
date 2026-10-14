include(OpenEyeSWIG)

set(OEMMPA_SWIG_EXPECTED_LIB_VARS
    OECHEM_LIBRARY OESYSTEM_LIBRARY OEPLATFORM_LIBRARY
    OEMATH_LIBRARY OEZSTD_LIBRARY
)

if(OEMMPA_HAS_OEMEDCHEM)
    list(APPEND OEMMPA_SWIG_EXPECTED_LIB_VARS OEMEDCHEM_LIBRARY)
endif()

set(OEMMPA_SWIG_FLAGS "-doxygen")
# SWIG %includes oedesalt's Desalter.h to wrap OEDESALT::Desalter into this
# module; the fetched headers are not on SWIG's default search path, so add
# them explicitly. (The generated wrapper's C++ compile already sees them
# transitively through the oemmpa target's link to oedesalt.)
list(APPEND OEMMPA_SWIG_FLAGS "-I${oedesalt_SOURCE_DIR}/include")
if(OEMMPA_HAS_DUCKDB)
    list(APPEND OEMMPA_SWIG_FLAGS "-DOEMMPA_HAS_DUCKDB=1")
else()
    list(APPEND OEMMPA_SWIG_FLAGS "-DOEMMPA_HAS_DUCKDB=0")
endif()

if(OEMMPA_HAS_OEMEDCHEM)
    list(APPEND OEMMPA_SWIG_FLAGS "-DOEMMPA_HAS_OEMEDCHEM=1")
else()
    list(APPEND OEMMPA_SWIG_FLAGS "-DOEMMPA_HAS_OEMEDCHEM=0")
endif()

openeye_add_swig_module(
    NAME oemmpa
    SWIG_FILE oemmpa.i
    LINK_LIBS oemmpa
    PYTHON_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/python/oemmpa
    STABLE_ABI ${OEMMPA_USE_STABLE_ABI}
    SWIG_FLAGS ${OEMMPA_SWIG_FLAGS}
    EXPECTED_LIB_VARS
        ${OEMMPA_SWIG_EXPECTED_LIB_VARS}
    INIT_PY ${CMAKE_SOURCE_DIR}/python/oemmpa/__init__.py
)

install(FILES
    ${CMAKE_SOURCE_DIR}/python/oemmpa/__main__.py
    ${CMAKE_SOURCE_DIR}/python/oemmpa/_analytics.py
    ${CMAKE_SOURCE_DIR}/python/oemmpa/_dataframe.py
    ${CMAKE_SOURCE_DIR}/python/oemmpa/_facade.py
    ${CMAKE_SOURCE_DIR}/python/oemmpa/_loading.py
    ${CMAKE_SOURCE_DIR}/python/oemmpa/_query.py
    ${CMAKE_SOURCE_DIR}/python/oemmpa/_results.py
    ${CMAKE_SOURCE_DIR}/python/oemmpa/_rgroup.py
    ${CMAKE_SOURCE_DIR}/python/oemmpa/_rule_environment.py
    ${CMAKE_SOURCE_DIR}/python/oemmpa/_series.py
    ${CMAKE_SOURCE_DIR}/python/oemmpa/_smiles_file.py
    ${CMAKE_SOURCE_DIR}/python/oemmpa/_storage.py
    ${CMAKE_SOURCE_DIR}/python/oemmpa/_transform.py
    ${CMAKE_SOURCE_DIR}/python/oemmpa/cli.py
    DESTINATION oemmpa
)

# Every build tree copies the extension into the one source-tree package
# directory, so the tree that built last owns what `import oemmpa` loads. A
# debug extension answers a statistics query a few ULP away from a release
# one, which is enough to make an acceptance digest disagree with itself.
#
# The helper stages the module from a POST_BUILD command, which runs only when
# the module actually relinks. That makes the documented recovery -- build the
# release tree last -- silently do nothing whenever the release tree is
# already up to date, leaving the other tree's extension in place. This target
# is a custom target, so it is always out of date and always restages; the
# copies are copy_if_different, so the cost when nothing moved is a stat.
#
# The marker's own name carries the configuration: a multi-config generator
# evaluates file(GENERATE) once per configuration, and CMake rejects one
# output path that would hold different content in each.
file(GENERATE
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/_build_type_$<CONFIG>.txt
    CONTENT "$<CONFIG>\n"
)
add_custom_target(oemmpa_python_stage ALL
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        $<TARGET_FILE:oemmpa_python>
        ${CMAKE_SOURCE_DIR}/python/oemmpa/$<TARGET_FILE_NAME:oemmpa_python>
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/_build_type_$<CONFIG>.txt
        ${CMAKE_SOURCE_DIR}/python/oemmpa/_build_type.txt
    DEPENDS oemmpa_python
    COMMENT "Staging the $<CONFIG> extension into python/oemmpa"
)
