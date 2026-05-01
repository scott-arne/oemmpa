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

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DuckDB
    REQUIRED_VARS DuckDB_INCLUDE_DIR DuckDB_LIBRARY
)

if(DuckDB_FOUND AND NOT TARGET DuckDB::DuckDB)
    add_library(DuckDB::DuckDB UNKNOWN IMPORTED)
    set_target_properties(DuckDB::DuckDB PROPERTIES
        IMPORTED_LOCATION "${DuckDB_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${DuckDB_INCLUDE_DIR}"
    )
endif()

mark_as_advanced(DuckDB_INCLUDE_DIR DuckDB_LIBRARY)
