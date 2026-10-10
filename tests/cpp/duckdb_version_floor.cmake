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

if(failures)
    list(JOIN failures "\n  " report)
    message(FATAL_ERROR "DuckDB version floor cases failed:\n  ${report}")
endif()
message(STATUS "DuckDB version floor: 6 cases passed (floor ${OEMMPA_DUCKDB_MIN_VERSION})")
