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

# Recorded in _build_info.py so the package can tell, at import, whether its
# DuckDB symbols could be coalesced with another copy's. A statically linked
# build exports none and is immune; a shared one is not, and the import guard
# only applies there (see python/oemmpa/_duckdb_runtime.py).
if(OEMMPA_HAS_DUCKDB AND OEMMPA_DUCKDB_STATIC)
    set(OEMMPA_DUCKDB_LINKAGE "STATIC")
elseif(OEMMPA_HAS_DUCKDB)
    set(OEMMPA_DUCKDB_LINKAGE "SHARED")
else()
    set(OEMMPA_DUCKDB_LINKAGE "NONE")
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
    EXTRA_BUILD_INFO "DUCKDB_LINKAGE = \"${OEMMPA_DUCKDB_LINKAGE}\""
)

# An editable build loads the extension from python/oemmpa, and the RPATH that
# openeye_add_swig_module() sets covers only that directory and the OpenEye
# runtime. Homebrew's libduckdb records its absolute install name, but the
# libduckdb.so in DuckDB's release archive has the bare SONAME libduckdb.so, so
# on Linux the import fails unless DuckDB's directory is on the RPATH. The entry
# goes on the target, not into the staged copy, because oemmpa_python_stage
# below re-copies the build-tree module over python/oemmpa on every build.
# Wheels are excluded: auditwheel and delocate vendor DuckDB and rewrite the
# RPATH, and a build machine's absolute path does not belong in a wheel.
if(NOT SKBUILD AND NOT WIN32 AND OEMMPA_DUCKDB_LINKAGE STREQUAL "SHARED")
    get_filename_component(_OEMMPA_DUCKDB_LIB_DIR "${DuckDB_LIBRARY}" DIRECTORY)
    set_property(TARGET oemmpa_python APPEND PROPERTY
        INSTALL_RPATH "${_OEMMPA_DUCKDB_LIB_DIR}")
endif()

# Cutting the extension's exports to its init symbol is the other half of the
# static link: an archive's symbols are still exported from the module that
# absorbs them, and an exported DuckDB symbol can be coalesced with the pip
# duckdb package's copy exactly as a dylib's can. Nothing dlopens anything out
# of this module, so one export is all it needs.
if(APPLE AND OEMMPA_HAS_DUCKDB AND OEMMPA_DUCKDB_STATIC)
    file(GENERATE
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/exported_symbols.txt
        CONTENT "_PyInit_$<TARGET_FILE_BASE_NAME:oemmpa_python>\n"
    )
    target_link_options(oemmpa_python PRIVATE
        -Wl,-exported_symbols_list,${CMAKE_CURRENT_BINARY_DIR}/exported_symbols.txt)
endif()

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
