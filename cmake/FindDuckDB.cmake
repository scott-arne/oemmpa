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

# DuckDB_VERSION comes from a DuckDBConfigVersion.cmake bound to the found
# library: the one beside it (Homebrew's layout), or one under a DUCKDB_ROOT
# hint that contains it (see oemmpa_duckdb_version_file_for). The official flat
# libduckdb bundle ships none, so DuckDB_VERSION stays empty and the caller
# treats the version as unknown. DuckDB_VERSION_FILE names the file the version
# was read from.
include("${CMAKE_CURRENT_LIST_DIR}/DuckDBVersion.cmake")
set(DuckDB_VERSION "")
set(DuckDB_VERSION_FILE "")
if(DuckDB_LIBRARY)
    oemmpa_duckdb_version_file_for("${DuckDB_LIBRARY}" DuckDB_VERSION_FILE
        "${DUCKDB_ROOT}" "$ENV{DUCKDB_ROOT}")
    if(NOT DuckDB_VERSION_FILE STREQUAL "")
        oemmpa_read_duckdb_version("${DuckDB_VERSION_FILE}" DuckDB_VERSION)
    endif()
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
