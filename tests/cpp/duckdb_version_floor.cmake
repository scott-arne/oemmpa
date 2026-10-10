# Exercises the DuckDB version parser and the configure-time floor on
# temporary DuckDBConfigVersion.cmake files.
cmake_minimum_required(VERSION 3.16)

include("${MODULE}")

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")

# The body Homebrew ships, so the parser is tested against the real shape,
# including the PACKAGE_VERSION_* lines it must not mistake for the version.
set(version_file_template [=[
set(PACKAGE_VERSION "@VERSION@")

if("${PACKAGE_VERSION}" VERSION_LESS "${PACKAGE_FIND_VERSION}")
    set(PACKAGE_VERSION_COMPATIBLE FALSE)
else()
    set(PACKAGE_VERSION_COMPATIBLE TRUE)
    if ("${PACKAGE_VERSION}" VERSION_EQUAL "${PACKAGE_FIND_VERSION}")
        set(PACKAGE_VERSION_EXACT TRUE)
    endif()
endif()
]=])

set(failures "")

# check_case(<name> <body variable, or empty for no file> <version> <status>)
#
# The body is passed by variable name so its ${...} text reaches the file
# verbatim.
function(check_case name body_var expected_version expected_status)
    set(path "${WORK_DIR}/${name}/DuckDBConfigVersion.cmake")
    if(NOT body_var STREQUAL "")
        file(WRITE "${path}" "${${body_var}}")
    endif()
    oemmpa_read_duckdb_version("${path}" actual_version)
    oemmpa_classify_duckdb_version("${actual_version}" actual_status)
    if(NOT actual_version STREQUAL expected_version OR NOT actual_status STREQUAL expected_status)
        list(APPEND failures "${name}: got '${actual_version}' (${actual_status}), expected '${expected_version}' (${expected_status})")
        set(failures "${failures}" PARENT_SCOPE)
    endif()
endfunction()

foreach(spec IN ITEMS "prefixed_below|v1.5.3|1.5.3|BELOW" "plain_below|1.5.3|1.5.3|BELOW"
                      "prefixed_above|v1.5.6|1.5.6|OK" "plain_at_floor|1.5.4|1.5.4|OK")
    string(REPLACE "|" ";" fields "${spec}")
    list(GET fields 0 name)
    list(GET fields 1 declared)
    list(GET fields 2 expected_version)
    list(GET fields 3 expected_status)
    string(REPLACE "@VERSION@" "${declared}" body "${version_file_template}")
    check_case("${name}" body "${expected_version}" "${expected_status}")
endforeach()

set(unversioned_body "set(PACKAGE_VERSION_COMPATIBLE TRUE)\n")
check_case(no_package_version unversioned_body "" UNKNOWN)
check_case(missing_file "" "" UNKNOWN)

# A version file only describes the library installed with it, so the lookup
# below is checked on temporary layouts: libraries are empty files, since only
# their paths matter.
set(layouts "${WORK_DIR}/layouts")

function(write_version_file directory declared)
    string(REPLACE "@VERSION@" "${declared}" body "${version_file_template}")
    file(WRITE "${directory}/cmake/DuckDB/DuckDBConfigVersion.cmake" "${body}")
endfunction()

# check_binding(<name> <library> <expected file or empty> <expected status> [<root>...])
function(check_binding name library expected_file expected_status)
    oemmpa_duckdb_version_file_for("${library}" actual_file ${ARGN})
    oemmpa_read_duckdb_version("${actual_file}" actual_version)
    oemmpa_classify_duckdb_version("${actual_version}" actual_status)
    if(NOT actual_file STREQUAL expected_file OR NOT actual_status STREQUAL expected_status)
        list(APPEND failures "${name}: got '${actual_file}' (${actual_status}), expected '${expected_file}' (${expected_status})")
        set(failures "${failures}" PARENT_SCOPE)
    endif()
endfunction()

# (a) The version file beside the library describes it.
file(WRITE "${layouts}/adjacent/lib/libduckdb.dylib" "")
write_version_file("${layouts}/adjacent/lib" v1.5.6)
check_binding(adjacent_file "${layouts}/adjacent/lib/libduckdb.dylib"
    "${layouts}/adjacent/lib/cmake/DuckDB/DuckDBConfigVersion.cmake" OK)

# (b) A root that contains the library describes it, here from lib64 so only
# the root rule can find the file.
file(WRITE "${layouts}/root/lib64/libduckdb.dylib" "")
write_version_file("${layouts}/root/lib" v1.5.6)
check_binding(containing_root "${layouts}/root/lib64/libduckdb.dylib"
    "${layouts}/root/lib/cmake/DuckDB/DuckDBConfigVersion.cmake" OK "${layouts}/root")

# The same root named through a symlink still contains the library.
file(CREATE_LINK "${layouts}/root" "${layouts}/root-alias" SYMBOLIC)
check_binding(containing_root_alias "${layouts}/root/lib64/libduckdb.dylib"
    "${layouts}/root-alias/lib/cmake/DuckDB/DuckDBConfigVersion.cmake" OK "${layouts}/root-alias")

# (c) and (d) A root that does not contain the library describes another
# DuckDB, whether its file would pass the floor or fail it.
file(WRITE "${layouts}/flat/libduckdb.dylib" "")
write_version_file("${layouts}/newer-root/lib" v1.5.6)
check_binding(unrelated_newer_root "${layouts}/flat/libduckdb.dylib" "" UNKNOWN "${layouts}/newer-root")
write_version_file("${layouts}/older-root/lib" v1.5.3)
check_binding(unrelated_older_root "${layouts}/flat/libduckdb.dylib" "" UNKNOWN "${layouts}/older-root")

# (e) A library reached through a symlink is described by the file beside the
# link (Homebrew's /opt/homebrew/lib), not only beside its target.
file(WRITE "${layouts}/cellar/1.5.6/lib/libduckdb.dylib" "")
file(MAKE_DIRECTORY "${layouts}/linked/lib")
file(CREATE_LINK "${layouts}/cellar/1.5.6/lib/libduckdb.dylib" "${layouts}/linked/lib/libduckdb.dylib" SYMBOLIC)
write_version_file("${layouts}/linked/lib" v1.5.6)
check_binding(symlinked_library "${layouts}/linked/lib/libduckdb.dylib"
    "${layouts}/linked/lib/cmake/DuckDB/DuckDBConfigVersion.cmake" OK)

if(failures)
    list(JOIN failures "\n  " report)
    message(FATAL_ERROR "DuckDB version floor cases failed:\n  ${report}")
endif()
message(STATUS "DuckDB version floor: 12 cases passed (floor ${OEMMPA_DUCKDB_MIN_VERSION})")
