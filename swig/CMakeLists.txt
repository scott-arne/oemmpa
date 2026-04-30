include(OpenEyeSWIG)

openeye_add_swig_module(
    NAME oemmpa
    SWIG_FILE oemmpa.i
    LINK_LIBS oemmpa
    PYTHON_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/python/oemmpa
    STABLE_ABI ${OEMMPA_USE_STABLE_ABI}
    SWIG_FLAGS "-doxygen"
    EXPECTED_LIB_VARS
        OECHEM_LIBRARY OESYSTEM_LIBRARY OEPLATFORM_LIBRARY
        OEMATH_LIBRARY OEZSTD_LIBRARY
    INIT_PY ${CMAKE_SOURCE_DIR}/python/oemmpa/__init__.py
)
