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

# oemmpa_duckdb_version_file_for(<library> <out_var> [<root>...])
#
# Set <out_var> to the DuckDBConfigVersion.cmake that declares the version of
# <library>, or to an empty string when none does. A version file describes
# only the DuckDB installed with it, so two places qualify: <library>'s own
# directory as given (<dir>/cmake/DuckDB/, which is Homebrew's layout read
# through its lib/ symlink), and <root>/lib/cmake/DuckDB/ for a <root> that
# contains <library> once both are resolved to real paths. A root that does
# not contain the library, such as a stale DUCKDB_ROOT beside a cached or
# explicitly supplied DuckDB_LIBRARY, describes some other DuckDB and is never
# read: its version could pass an old library or fail a valid one. Files that
# declare no version are skipped.
function(oemmpa_duckdb_version_file_for library out_var)
    get_filename_component(library_dir "${library}" DIRECTORY)
    set(candidates "${library_dir}/cmake/DuckDB/DuckDBConfigVersion.cmake")
    get_filename_component(library_real "${library}" REALPATH)
    foreach(root IN LISTS ARGN)
        if(root STREQUAL "")
            continue()
        endif()
        get_filename_component(root_real "${root}" REALPATH)
        if(root_real STREQUAL "/")
            set(root_prefix "/")
        else()
            set(root_prefix "${root_real}/")
        endif()
        string(FIND "${library_real}" "${root_prefix}" position)
        if(position EQUAL 0)
            list(APPEND candidates "${root}/lib/cmake/DuckDB/DuckDBConfigVersion.cmake")
        endif()
    endforeach()
    set(found "")
    foreach(candidate IN LISTS candidates)
        oemmpa_read_duckdb_version("${candidate}" version)
        if(NOT version STREQUAL "")
            set(found "${candidate}")
            break()
        endif()
    endforeach()
    set(${out_var} "${found}" PARENT_SCOPE)
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
