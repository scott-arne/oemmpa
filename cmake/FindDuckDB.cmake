find_path(DuckDB_INCLUDE_DIR
    NAMES duckdb.hpp
    HINTS
        ${DUCKDB_ROOT}
        $ENV{DUCKDB_ROOT}
    PATHS
        /opt/homebrew/opt/duckdb
        /opt/homebrew
        /usr/local/opt/duckdb
        /usr/local
    PATH_SUFFIXES include
)

find_library(DuckDB_LIBRARY
    NAMES duckdb
    HINTS
        ${DUCKDB_ROOT}
        $ENV{DUCKDB_ROOT}
    PATHS
        /opt/homebrew/opt/duckdb
        /opt/homebrew
        /usr/local/opt/duckdb
        /usr/local
    PATH_SUFFIXES lib
)

# DuckDB_VERSION comes from the DuckDBConfigVersion.cmake beside the found
# library (Homebrew's layout), else from one under the DUCKDB_ROOT hints. The
# official flat libduckdb bundle ships none, so DuckDB_VERSION stays empty and
# the caller treats the version as unknown. DuckDB_VERSION_FILE names the file
# the version was read from.
include("${CMAKE_CURRENT_LIST_DIR}/DuckDBVersion.cmake")
set(DuckDB_VERSION "")
set(DuckDB_VERSION_FILE "")
if(DuckDB_LIBRARY)
    get_filename_component(_duckdb_library_dir "${DuckDB_LIBRARY}" DIRECTORY)
    set(_duckdb_version_candidates
        "${_duckdb_library_dir}/cmake/DuckDB/DuckDBConfigVersion.cmake")
    foreach(_duckdb_root IN ITEMS "${DUCKDB_ROOT}" "$ENV{DUCKDB_ROOT}")
        if(NOT _duckdb_root STREQUAL "")
            list(APPEND _duckdb_version_candidates
                "${_duckdb_root}/lib/cmake/DuckDB/DuckDBConfigVersion.cmake")
        endif()
    endforeach()
    foreach(_duckdb_candidate IN LISTS _duckdb_version_candidates)
        oemmpa_read_duckdb_version("${_duckdb_candidate}" _duckdb_version)
        if(NOT _duckdb_version STREQUAL "")
            set(DuckDB_VERSION "${_duckdb_version}")
            set(DuckDB_VERSION_FILE "${_duckdb_candidate}")
            break()
        endif()
    endforeach()
    unset(_duckdb_library_dir)
    unset(_duckdb_version_candidates)
    unset(_duckdb_root)
    unset(_duckdb_candidate)
    unset(_duckdb_version)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DuckDB
    REQUIRED_VARS DuckDB_INCLUDE_DIR DuckDB_LIBRARY
    VERSION_VAR DuckDB_VERSION
)

if(DuckDB_FOUND AND NOT TARGET DuckDB::DuckDB)
    add_library(DuckDB::DuckDB UNKNOWN IMPORTED)
    set_target_properties(DuckDB::DuckDB PROPERTIES
        IMPORTED_LOCATION "${DuckDB_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${DuckDB_INCLUDE_DIR}"
    )
endif()

mark_as_advanced(DuckDB_INCLUDE_DIR DuckDB_LIBRARY)
