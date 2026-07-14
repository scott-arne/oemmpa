cmake_minimum_required(VERSION 3.16)

if(POLICY CMP0144)
    cmake_policy(SET CMP0144 NEW)
endif()

project(oemmpa VERSION 2.2.2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(OEMMPA_BUILD_TESTS "Build tests" ON)
option(OEMMPA_BUILD_PYTHON "Build Python bindings" ON)
option(OEMMPA_UNIVERSAL2 "Build universal2 binary for macOS" OFF)
option(OEMMPA_USE_STABLE_ABI "Use Python stable ABI" ON)
option(OEMMPA_BUILD_DUCKDB "Build optional DuckDB-backed storage when DuckDB is available" OFF)

if(APPLE AND OEMMPA_UNIVERSAL2)
    set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64" CACHE STRING "Build architectures" FORCE)
elseif(APPLE)
    unset(CMAKE_OSX_ARCHITECTURES CACHE)
endif()

if(APPLE AND NOT CMAKE_OSX_DEPLOYMENT_TARGET)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0" CACHE STRING "Minimum macOS version")
endif()

# Dependencies
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# The cmake-openeye modules (FindOpenEye, FindOpenEyePython, OpenEyeSWIG) are
# vendored under cmake/cmake-openeye so configure never depends on live GitHub
# access. See cmake/cmake-openeye/README.md for the pinned version and refresh
# instructions.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake-openeye")

# Desalting is provided by the standalone oedesalt library, fetched and linked
# statically. Its own build reuses our OpenEye/Python discovery (the imported
# targets are created idempotently). The fetch needs GitHub access at configure
# time; behind the corporate proxy, export HTTP(S)_PROXY before configuring
# (mirrors how oedesalt itself fetches cmake-openeye). We disable oedesalt's
# tests and Python bindings — oemmpa links only the C++ static library and wraps
# the Desalter through its own SWIG module.
include(FetchContent)
set(OEDESALT_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(OEDESALT_BUILD_PYTHON OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
    oedesalt
    GIT_REPOSITORY https://github.com/scott-arne/oedesalt.git
    GIT_TAG v0.1.0
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(oedesalt)

if(OEMMPA_BUILD_PYTHON)
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
    # Python bindings cannot be built without SWIG, so fail fast with an
    # explicit message rather than silently producing a build with no _oemmpa
    # extension. Opt out with -DOEMMPA_BUILD_PYTHON=OFF for a C++-only build.
    find_package(SWIG 4.0)
    if(NOT SWIG_FOUND)
        message(FATAL_ERROR
            "OEMMPA_BUILD_PYTHON is ON but SWIG (>= 4.0) was not found. "
            "Install SWIG to build the Python bindings, or configure with "
            "-DOEMMPA_BUILD_PYTHON=OFF for a C++-only build.")
    endif()
endif()

if(OPENEYE_USE_SHARED)
    find_package(OpenEyePython REQUIRED)
endif()

find_package(OpenEye REQUIRED)

if(OEMMPA_BUILD_DUCKDB)
    find_package(DuckDB QUIET)
endif()

if(OEMMPA_BUILD_DUCKDB AND DuckDB_FOUND)
    set(OEMMPA_HAS_DUCKDB TRUE)
else()
    set(OEMMPA_HAS_DUCKDB FALSE)
    if(OEMMPA_BUILD_DUCKDB)
        message(STATUS "DuckDB not found; DuckDB-backed storage is disabled")
    endif()
endif()

if(TARGET OpenEye::OEMedChem)
    set(OEMMPA_HAS_OEMEDCHEM TRUE)
else()
    set(OEMMPA_HAS_OEMEDCHEM FALSE)
endif()

# Library
set(OEMMPA_SOURCES
    src/oemmpa.cpp
    src/DatabaseSummary.cpp
    src/RuleEnvironmentStatistics.cpp
    src/MoleculeRecord.cpp
    src/LoadReport.cpp
    src/Fragmentation.cpp
    src/FragmentationStrategy.cpp
    src/Fragmenter.cpp
    src/EnvironmentFingerprint.cpp
    src/QueryEnvironment.cpp
    src/FragmentationMethod.cpp
    src/MCSCommon.cpp
    src/DMCSSMethod.cpp
    src/Analyzer.cpp
    src/MatchedPair.cpp
    src/Transform.cpp
    src/TransformApplication.cpp
    src/QueryOptions.cpp
    src/PairScoring.cpp
    src/MemoryIndex.cpp
    src/ThreadSupport.cpp
)

if(OEMMPA_HAS_DUCKDB)
    list(APPEND OEMMPA_SOURCES src/DuckDBStore.cpp)
endif()

if(OEMMPA_HAS_OEMEDCHEM)
    list(APPEND OEMMPA_SOURCES src/OEMedChemMethod.cpp)
endif()

add_library(oemmpa ${OEMMPA_SOURCES})

set_target_properties(oemmpa PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(oemmpa
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(oemmpa
    PUBLIC
        OpenEye::OEChem
        # oedesalt is a fetched static library. Scope the link to the build tree
        # with BUILD_INTERFACE so it (and its own OpenEye usage requirements)
        # applies to the in-tree SWIG module and tests, while staying out of the
        # installed oemmpa-targets export set — oedesalt manages its own export.
        $<BUILD_INTERFACE:oedesalt>
)

if(OEMMPA_HAS_OEMEDCHEM)
    target_compile_definitions(oemmpa PUBLIC OEMMPA_HAS_OEMEDCHEM=1)
    target_link_libraries(oemmpa PUBLIC OpenEye::OEMedChem)
else()
    target_compile_definitions(oemmpa PUBLIC OEMMPA_HAS_OEMEDCHEM=0)
endif()

if(OEMMPA_HAS_DUCKDB)
    target_compile_definitions(oemmpa PUBLIC OEMMPA_HAS_DUCKDB=1)
    target_link_libraries(oemmpa PUBLIC DuckDB::DuckDB)
else()
    target_compile_definitions(oemmpa PUBLIC OEMMPA_HAS_DUCKDB=0)
endif()

# Tests
if(OEMMPA_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Python bindings. SWIG is REQUIRED above when OEMMPA_BUILD_PYTHON is ON, so
# reaching here guarantees it was found.
if(OEMMPA_BUILD_PYTHON)
    add_subdirectory(swig)
endif()

# Install
install(TARGETS oemmpa
    EXPORT oemmpa-targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(DIRECTORY include/oemmpa DESTINATION include)

install(EXPORT oemmpa-targets
    FILE OEMMPATargets.cmake
    NAMESPACE OEMMPA::
    DESTINATION lib/cmake/oemmpa
)
