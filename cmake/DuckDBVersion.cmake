# Configure-time DuckDB version floor.
#
# DuckDB's CMake package declares its version with a leading "v" (Homebrew
# writes set(PACKAGE_VERSION "v1.5.6")), which CMake's VERSION_* comparisons
# cannot order. The version is therefore read from the file's text with the
# "v" stripped, and the floor is checked by the caller after find_package
# rather than requested from it.

include_guard(GLOBAL)

# The oldest DuckDB the storage layer is written against (sequence handling
# and transaction rules in src/DuckDBStore.cpp assume it).
set(OEMMPA_DUCKDB_MIN_VERSION "1.5.4")

# oemmpa_read_duckdb_version(<path> <out_var>)
#
# Set <out_var> to the version <path>, a DuckDBConfigVersion.cmake, assigns to
# PACKAGE_VERSION, without a leading "v". <out_var> is empty when the file is
# missing or assigns no literal version. The file is parsed rather than
# included so that none of its code runs in the caller's scope.
function(oemmpa_read_duckdb_version path out_var)
    set(version "")
    if(EXISTS "${path}" AND NOT IS_DIRECTORY "${path}")
        file(STRINGS "${path}" lines REGEX "PACKAGE_VERSION[ \t]")
        foreach(line IN LISTS lines)
            if(line MATCHES "^[ \t]*[Ss][Ee][Tt][ \t]*\\([ \t]*PACKAGE_VERSION[ \t]+\"?[vV]?([0-9][0-9A-Za-z.+-]*)\"?[ \t]*\\)")
                set(version "${CMAKE_MATCH_1}")
                break()
            endif()
        endforeach()
    endif()
    set(${out_var} "${version}" PARENT_SCOPE)
endfunction()

# oemmpa_classify_duckdb_version(<version> <out_var>)
#
# Set <out_var> to UNKNOWN when <version> is empty, BELOW when it is older
# than OEMMPA_DUCKDB_MIN_VERSION, and OK otherwise.
function(oemmpa_classify_duckdb_version version out_var)
    if(version STREQUAL "")
        set(status UNKNOWN)
    elseif(version VERSION_LESS OEMMPA_DUCKDB_MIN_VERSION)
        set(status BELOW)
    else()
        set(status OK)
    endif()
    set(${out_var} "${status}" PARENT_SCOPE)
endfunction()
