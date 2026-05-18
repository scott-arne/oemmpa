include(OpenEyeSWIG)

set(OEMMPA_SWIG_EXPECTED_LIB_VARS
    OECHEM_LIBRARY OESYSTEM_LIBRARY OEPLATFORM_LIBRARY
    OEMATH_LIBRARY OEZSTD_LIBRARY
)

if(OEMMPA_HAS_OEMEDCHEM)
    list(APPEND OEMMPA_SWIG_EXPECTED_LIB_VARS OEMEDCHEM_LIBRARY)
endif()

set(OEMMPA_SWIG_FLAGS "-doxygen")
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
    ${CMAKE_SOURCE_DIR}/python/oemmpa/_smiles_file.py
    ${CMAKE_SOURCE_DIR}/python/oemmpa/_storage.py
    ${CMAKE_SOURCE_DIR}/python/oemmpa/_transform.py
    ${CMAKE_SOURCE_DIR}/python/oemmpa/cli.py
    DESTINATION oemmpa
)
