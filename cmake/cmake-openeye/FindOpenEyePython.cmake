# FindOpenEyePython.cmake
# Auto-discover OpenEye shared libraries from the openeye-toolkits Python package
#
# This module queries the installed openeye-toolkits Python package to find
# the shared library directory and marketing version. It sets variables that
# FindOpenEye.cmake and OpenEyeSWIG.cmake consume
# (OPENEYE_RUNTIME_LIB_DIR for RPATH, OpenEyePython_PLATFORM).
#
# Requirements:
#   - Python3_EXECUTABLE must be set (call find_package(Python3) first)
#   - openeye-toolkits must be installed in the Python environment
#
# The following variables are set:
#   OpenEyePython_FOUND            - TRUE if openeye-toolkits was found
#   OpenEyePython_LIB_DIR          - Absolute path to the shared library directory
#   OpenEyePython_VERSION          - Marketing version (e.g., "2025.2.1")
#   OpenEyePython_PLATFORM         - Platform subdirectory name (e.g., "osx-arm64-14-clang-15.0")
#   OPENEYE_LIB_DIR                - Set on POSIX for FindOpenEye.cmake shared-link discovery
#   OPENEYE_RUNTIME_LIB_DIR        - Set for OpenEyeSWIG.cmake runtime/RPATH handling
#   OPENEYE_USE_SHARED             - Set ON on POSIX for openeye-toolkits shared-library builds
#   OPENEYE_TOOLKITS_VERSION       - Set for _build_info.py generation

# Check that Python3 is available
if(NOT Python3_EXECUTABLE)
    message(FATAL_ERROR "FindOpenEyePython requires Python3_EXECUTABLE to be set. "
        "Call find_package(Python3 COMPONENTS Interpreter) before find_package(OpenEyePython).")
endif()

# Query the library directory and platform from openeye-toolkits.
# On Windows, DLLs live flat in openeye/libs/ with no arch subdirectory,
# and upstream libs.FindOpenEyeDLLSDirectory() is broken (iterates characters
# of a string), so probe openeye.libs.__file__ directly on that platform.
if(WIN32)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c
            "import openeye.libs, os; d = os.path.dirname(openeye.libs.__file__); print(d); print(os.path.basename(d))"
        OUTPUT_VARIABLE _OE_PYTHON_OUTPUT
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE _OE_PYTHON_RESULT
    )
else()
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c
            "from openeye import libs; import os; d = libs.FindOpenEyeDLLSDirectory(); print(d); print(os.path.basename(d))"
        OUTPUT_VARIABLE _OE_PYTHON_OUTPUT
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE _OE_PYTHON_RESULT
    )
endif()

if(NOT _OE_PYTHON_RESULT EQUAL 0)
    message(STATUS "FindOpenEyePython: openeye-toolkits Python package not found or FindOpenEyeDLLSDirectory() failed")
    set(OpenEyePython_FOUND FALSE)
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(OpenEyePython
        REQUIRED_VARS OpenEyePython_LIB_DIR
        FAIL_MESSAGE "openeye-toolkits Python package not found. Install with: pip install openeye-toolkits"
    )
    return()
endif()

# Parse the two-line output: line 1 = full lib dir, line 2 = platform basename
string(REPLACE "\n" ";" _OE_PYTHON_LINES "${_OE_PYTHON_OUTPUT}")
list(GET _OE_PYTHON_LINES 0 _OE_LIB_DIR)
list(GET _OE_PYTHON_LINES 1 _OE_PLATFORM)

message(STATUS "FindOpenEyePython: Library directory: ${_OE_LIB_DIR}")
message(STATUS "FindOpenEyePython: Platform: ${_OE_PLATFORM}")

# Query the marketing version
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c
        "from openeye import oechem; print(oechem.OEToolkitsGetRelease())"
    OUTPUT_VARIABLE _OE_MARKETING_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    RESULT_VARIABLE _OE_VERSION_RESULT
)

if(NOT _OE_VERSION_RESULT EQUAL 0)
    message(STATUS "FindOpenEyePython: Could not query OEToolkitsGetRelease(), version will be empty")
    set(_OE_MARKETING_VERSION "")
endif()

if(_OE_MARKETING_VERSION)
    message(STATUS "FindOpenEyePython: Marketing version: ${_OE_MARKETING_VERSION}")
endif()

# Validate that the library directory contains shared libraries
if(APPLE)
    file(GLOB _OE_SHARED_LIBS "${_OE_LIB_DIR}/*.dylib")
elseif(WIN32)
    file(GLOB _OE_SHARED_LIBS "${_OE_LIB_DIR}/*.dll")
else()
    file(GLOB _OE_SHARED_LIBS "${_OE_LIB_DIR}/*.so")
endif()

if(NOT _OE_SHARED_LIBS)
    message(STATUS "FindOpenEyePython: No shared libraries found in ${_OE_LIB_DIR}")
    set(OpenEyePython_FOUND FALSE)
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(OpenEyePython
        REQUIRED_VARS OpenEyePython_LIB_DIR
        FAIL_MESSAGE "OpenEye shared libraries not found in ${_OE_LIB_DIR}"
    )
    return()
endif()

# Set output variables
set(OpenEyePython_LIB_DIR "${_OE_LIB_DIR}" CACHE PATH "OpenEye Python library directory")
set(OpenEyePython_VERSION "${_OE_MARKETING_VERSION}")
set(OpenEyePython_PLATFORM "${_OE_PLATFORM}")

# Set variables for FindOpenEye.cmake / OpenEyeSWIG.cmake consumption.
# OPENEYE_RUNTIME_LIB_DIR points at the wheel's openeye/libs/ directory (contains
# versioned .so/.dylib on POSIX, flat .dll on Windows). It is consumed by
# OpenEyeSWIG.cmake to compute RPATH on POSIX; Windows resolves DLLs at Python
# import time via openeye.libs' os.add_dll_directory() side effect.
set(OPENEYE_RUNTIME_LIB_DIR "${_OE_LIB_DIR}" CACHE PATH "Runtime shared-library directory from openeye-toolkits Python package" FORCE)
set(OPENEYE_TOOLKITS_VERSION "${_OE_MARKETING_VERSION}" CACHE STRING "OpenEye toolkits marketing version from Python package")
if(NOT WIN32)
    # POSIX openeye-toolkits ships linkable versioned .so/.dylib files in the
    # runtime library directory, so it is also the shared link-time directory.
    set(OPENEYE_LIB_DIR "${_OE_LIB_DIR}" CACHE PATH "OpenEye shared-library directory from openeye-toolkits Python package" FORCE)
    set(OPENEYE_USE_SHARED ON CACHE BOOL "Prefer shared OpenEye libraries from openeye-toolkits" FORCE)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenEyePython
    REQUIRED_VARS
        OpenEyePython_LIB_DIR
    VERSION_VAR
        OpenEyePython_VERSION
    FAIL_MESSAGE
        "openeye-toolkits Python package not found. Install with: pip install openeye-toolkits"
)

mark_as_advanced(
    OpenEyePython_LIB_DIR
)
