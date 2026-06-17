# FindOpenEye.cmake
# Find OpenEye Toolkits installation
#
# This module finds the OpenEye C++ toolkits and creates imported targets.
#
# User can set these environment variables or CMake variables:
#   OPENEYE_ROOT or OE_DIR - Root directory of OpenEye installation (for headers)
#   OPENEYE_LIB_DIR        - Override link-time library directory
#   OPENEYE_RUNTIME_LIB_DIR - Runtime shared-library directory from openeye-toolkits
#
# Options:
#   OPENEYE_USE_SHARED - Prefer shared libraries over static (default: OFF)
#                        Set to ON for wheels that depend on openeye-toolkits
#
# The following imported targets are created:
#   OpenEye::OEChem     - OEChem library (main chemistry library)
#   OpenEye::OESystem   - OESystem library
#   OpenEye::OEPlatform - OEPlatform library
#   OpenEye::OEMath     - OEMath library
#   OpenEye::OEGraphSim - OEGraphSim library (if available)
#   OpenEye::OEMedChem  - OEMedChem library (if available)
#   OpenEye::OEBio      - OEBio library (if available)
#   OpenEye::OEGrid     - OEGrid library (if available)
#   OpenEye::OEFizzChem - OEFizzChem library (if available)
#   OpenEye::zstd       - Bundled zstd library (if available)
#   OpenEye::OEOpt          - OEOpt library (if available)
#   OpenEye::OEMolPotential - OEMolPotential library (if available)
#   OpenEye::OEHermite      - OEHermite library (if available)
#   OpenEye::OEShape        - OEShape library (if available)
#   OpenEye::OEZap          - OEZap library (if available)
#   OpenEye::OESpicoli      - OESpicoli library (if available)
#   OpenEye::OESiteHopper   - OESiteHopper library (if available)
#   OpenEye::OEMMFF         - OEMMFF library (if available)
#   OpenEye::OEFF           - OEFF library (if available)
#   OpenEye::OESmirnoff     - OESmirnoff library (if available)
#   OpenEye::OEAmber        - OEAmber library (if available)
#   OpenEye::OEAM1          - OEAM1 library (if available)
#   OpenEye::OEAM1BCC       - OEAM1BCC library (if available)
#   OpenEye::OESzybki       - OESzybki library (if available)
#   OpenEye::OEQuacpac      - OEQuacpac library (if available)
#   OpenEye::OEOmega2       - OEOmega2 library (if available)
#   OpenEye::OESheffield    - OESheffield library (if available)
#   OpenEye::OESpruce       - OESpruce library (if available)
#   OpenEye::OEDepict       - OEDepict library (if available)
#   OpenEye::OEIUPAC        - OEIUPAC library (if available)
#
# The following variables are set:
#   OpenEye_FOUND              - TRUE if OpenEye was found
#   OpenEye_VERSION            - Version string (e.g., "4.3.0.1")
#   OpenEye_LIBRARY_TYPE       - SHARED or STATIC
#   OpenEye_GraphSim_FOUND     - TRUE if OEGraphSim was found
#   OpenEye_MedChem_FOUND      - TRUE if OEMedChem was found
#   OpenEye_Bio_FOUND          - TRUE if OEBio was found
#   OpenEye_Grid_FOUND         - TRUE if OEGrid was found
#   OpenEye_Opt_FOUND          - TRUE if OEOpt was found
#   OpenEye_MolPotential_FOUND - TRUE if OEMolPotential was found
#   OpenEye_Hermite_FOUND      - TRUE if OEHermite was found
#   OpenEye_Shape_FOUND        - TRUE if OEShape was found
#   OpenEye_Zap_FOUND          - TRUE if OEZap was found
#   OpenEye_Spicoli_FOUND      - TRUE if OESpicoli was found
#   OpenEye_SiteHopper_FOUND   - TRUE if OESiteHopper was found
#   OpenEye_MMFF_FOUND         - TRUE if OEMMFF was found
#   OpenEye_FF_FOUND           - TRUE if OEFF was found
#   OpenEye_Szybki_FOUND       - TRUE if OESzybki was found
#   OpenEye_Quacpac_FOUND      - TRUE if OEQuacpac was found
#   OpenEye_Omega2_FOUND       - TRUE if OEOmega2 was found
#   OpenEye_Sheffield_FOUND    - TRUE if OESheffield was found
#   OpenEye_Spruce_FOUND       - TRUE if OESpruce was found
#   OpenEye_Depict_FOUND       - TRUE if OEDepict was found
#   OpenEye_IUPAC_FOUND        - TRUE if OEIUPAC was found

option(OPENEYE_USE_SHARED "Prefer shared OpenEye libraries for dynamic linking" OFF)
set(OPENEYE_LIB_DIR "" CACHE PATH "Override OpenEye library directory for link-time discovery (e.g., SDK lib/)")
set(OPENEYE_RUNTIME_LIB_DIR "" CACHE PATH "Runtime shared-library directory (e.g., openeye-toolkits wheel openeye/libs/); used for RPATH on POSIX")

if(OPENEYE_RUNTIME_LIB_DIR)
    message(STATUS "OpenEye: Runtime library directory: ${OPENEYE_RUNTIME_LIB_DIR}")
endif()

# Windows ships only static .libs in the SDK, and the openeye-toolkits Python
# wheel has no MSVC import libraries at all — shared linking is not possible
# on Windows. Force OPENEYE_USE_SHARED OFF so the rest of this module can
# assume "shared mode" implies POSIX.
if(WIN32 AND OPENEYE_USE_SHARED)
    message(STATUS "OpenEye: OPENEYE_USE_SHARED is not supported on Windows (SDK is static-only); forcing OFF")
    set(OPENEYE_USE_SHARED OFF CACHE BOOL "" FORCE)
endif()

# Warn when using shared mode without a library directory override. POSIX uses
# the wheel's versioned libraries via OPENEYE_LIB_DIR or OPENEYE_RUNTIME_LIB_DIR.
if(OPENEYE_USE_SHARED AND NOT OPENEYE_LIB_DIR AND NOT OPENEYE_RUNTIME_LIB_DIR)
    message(WARNING "OPENEYE_USE_SHARED is ON but neither OPENEYE_LIB_DIR nor OPENEYE_RUNTIME_LIB_DIR is set. "
        "Shared library discovery may fail without an explicit library directory. "
        "Consider using FindOpenEyePython.cmake to auto-detect the library directory.")
endif()

# Look for the include directory. System fallbacks (/opt/openeye, /usr/local)
# are consulted only when no explicit hint was given, for the same reason as
# _OPENEYE_SYSTEM_LIB_FALLBACKS below: a host with a stray full SDK at
# /opt/openeye must not be able to override a caller who pointed us elsewhere.
if(OPENEYE_ROOT OR OE_DIR OR DEFINED ENV{OPENEYE_ROOT} OR DEFINED ENV{OE_DIR})
    set(_OPENEYE_SYSTEM_INCLUDE_FALLBACKS)
else()
    set(_OPENEYE_SYSTEM_INCLUDE_FALLBACKS /opt/openeye/include /usr/local/openeye/include)
endif()
find_path(OPENEYE_INCLUDE_DIR
    NAMES openeye.h
    PATHS
        ${OPENEYE_ROOT}/include
        ${OE_DIR}/include
        $ENV{OPENEYE_ROOT}/include
        $ENV{OE_DIR}/include
        ${_OPENEYE_SYSTEM_INCLUDE_FALLBACKS}
    PATH_SUFFIXES openeye
)

# Get the library directory - use override if provided, otherwise derive from include path
if(OPENEYE_LIB_DIR)
    message(STATUS "OpenEye: Using library directory override: ${OPENEYE_LIB_DIR}")
    set(_OPENEYE_LIB_SEARCH_PATHS ${OPENEYE_LIB_DIR})
elseif(OPENEYE_USE_SHARED AND OPENEYE_RUNTIME_LIB_DIR AND NOT WIN32)
    message(STATUS "OpenEye: Using runtime library directory for POSIX shared-link discovery: ${OPENEYE_RUNTIME_LIB_DIR}")
    set(_OPENEYE_LIB_SEARCH_PATHS ${OPENEYE_RUNTIME_LIB_DIR})
elseif(OPENEYE_INCLUDE_DIR)
    get_filename_component(_DEFAULT_LIB_DIR "${OPENEYE_INCLUDE_DIR}/../lib" ABSOLUTE)
    set(_OPENEYE_LIB_SEARCH_PATHS ${_DEFAULT_LIB_DIR})
endif()

# Set library search order based on preference (save/restore to not affect other finds)
set(_SAVED_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
if(OPENEYE_USE_SHARED)
    # POSIX-only (Windows coerced OFF above). Prefer dylib/so over static.
    if(APPLE)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .dylib .a)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES .so .a)
    endif()
    message(STATUS "OpenEye: Preferring shared libraries for dynamic linking")
endif()

# System fallback dirs are consulted only when no explicit hint was given.
# Otherwise a caller pointing at a partial/overlay SDK via OPENEYE_ROOT could
# silently pick up stray libraries from /opt/openeye/lib and mix-and-match
# the two trees — which broke test_optional_target_closure on CI runners
# that happened to have a full SDK installed at /opt/openeye.
if(OPENEYE_ROOT OR OE_DIR OR OPENEYE_LIB_DIR OR OPENEYE_RUNTIME_LIB_DIR
        OR DEFINED ENV{OPENEYE_ROOT} OR DEFINED ENV{OE_DIR})
    set(_OPENEYE_SYSTEM_LIB_FALLBACKS)
else()
    set(_OPENEYE_SYSTEM_LIB_FALLBACKS /opt/openeye/lib /usr/local/openeye/lib)
endif()

# Helper macro to find OpenEye library, handling versioned names (e.g., liboechem-4.3.0.1.dylib)
macro(find_openeye_library VAR_NAME LIB_NAME)
    # First try to find versioned shared library in the override directory
    # (openeye-toolkits Python package). POSIX-only; OPENEYE_USE_SHARED is
    # always OFF on Windows.
    if(OPENEYE_USE_SHARED)
        if(OPENEYE_LIB_DIR)
            set(_OPENEYE_VERSIONED_SEARCH_DIR "${OPENEYE_LIB_DIR}")
        elseif(OPENEYE_RUNTIME_LIB_DIR AND NOT WIN32)
            set(_OPENEYE_VERSIONED_SEARCH_DIR "${OPENEYE_RUNTIME_LIB_DIR}")
        else()
            set(_OPENEYE_VERSIONED_SEARCH_DIR "")
        endif()
    endif()
    if(_OPENEYE_VERSIONED_SEARCH_DIR)
        if(APPLE)
            file(GLOB _VERSIONED_LIB "${_OPENEYE_VERSIONED_SEARCH_DIR}/lib${LIB_NAME}-*.dylib")
        else()
            file(GLOB _VERSIONED_LIB "${_OPENEYE_VERSIONED_SEARCH_DIR}/lib${LIB_NAME}-*.so")
        endif()
        if(_VERSIONED_LIB)
            # Get the first match (should only be one)
            list(GET _VERSIONED_LIB 0 ${VAR_NAME})
            message(STATUS "OpenEye: Found versioned ${LIB_NAME}: ${${VAR_NAME}}")
        endif()
    endif()

    # Fall back to standard find_library if versioned library not found
    if(NOT ${VAR_NAME})
        find_library(${VAR_NAME}
            NAMES ${LIB_NAME}
            PATHS
                ${_OPENEYE_LIB_SEARCH_PATHS}
                ${OPENEYE_ROOT}/lib
                ${OE_DIR}/lib
                $ENV{OPENEYE_ROOT}/lib
                $ENV{OE_DIR}/lib
                ${_OPENEYE_SYSTEM_LIB_FALLBACKS}
            NO_DEFAULT_PATH
        )
    endif()
endmacro()

# Find required libraries
find_openeye_library(OECHEM_LIBRARY oechem)
find_openeye_library(OESYSTEM_LIBRARY oesystem)
find_openeye_library(OEPLATFORM_LIBRARY oeplatform)
find_openeye_library(OEMATH_LIBRARY oemath)

# Find optional libraries
find_openeye_library(OEGRAPHSIM_LIBRARY oegraphsim)
find_openeye_library(OEMEDCHEM_LIBRARY oemedchem)
find_openeye_library(OEBIO_LIBRARY oebio)
find_openeye_library(OEGRID_LIBRARY oegrid)
find_openeye_library(OEFIZZCHEM_LIBRARY oefizzchem)

# v1.1.0: additional library discovery for geometry/optimization
find_openeye_library(OEOPT_LIBRARY oeopt)
find_openeye_library(OEMOLPOTENTIAL_LIBRARY oemolpotential)
find_openeye_library(OEHERMITE_LIBRARY oehermite)
find_openeye_library(OESHAPE_LIBRARY oeshape)
find_openeye_library(OEZAP_LIBRARY oezap)
find_openeye_library(OESPICOLI_LIBRARY oespicoli)
find_openeye_library(OESITEHOPPER_LIBRARY oesitehopper)
find_openeye_library(OEMMFF_LIBRARY oemmff)
find_openeye_library(OEFF_LIBRARY oeff)
find_openeye_library(OESMIRNOFF_LIBRARY oesmirnoff)
find_openeye_library(OEAMBER_LIBRARY oeamber)
find_openeye_library(OEAM1_LIBRARY oeam1)
find_openeye_library(OEAM1BCC_LIBRARY oeam1bcc)
find_openeye_library(OESZYBKI_LIBRARY oeszybki)
find_openeye_library(OEQUACPAC_LIBRARY oequacpac)
find_openeye_library(OEOMEGA2_LIBRARY oeomega2)
find_openeye_library(OESHEFFIELD_LIBRARY oesheffield)
find_openeye_library(OESPRUCE_LIBRARY oespruce)
find_openeye_library(OEDEPICT_LIBRARY oedepict)
find_openeye_library(OEIUPAC_LIBRARY oeiupac)

# Find bundled zstd library (OpenEye bundles this) - uses different naming.
# Versioned-glob is POSIX-only (OPENEYE_USE_SHARED is always OFF on Windows).
# On Windows the SDK ships zstd_static.lib; find_library below picks it up so
# oeplatform's gzstd_* symbols resolve.
if(OPENEYE_USE_SHARED)
    if(OPENEYE_LIB_DIR)
        set(_OPENEYE_ZSTD_SEARCH_DIR "${OPENEYE_LIB_DIR}")
    elseif(OPENEYE_RUNTIME_LIB_DIR AND NOT WIN32)
        set(_OPENEYE_ZSTD_SEARCH_DIR "${OPENEYE_RUNTIME_LIB_DIR}")
    else()
        set(_OPENEYE_ZSTD_SEARCH_DIR "")
    endif()
endif()
if(_OPENEYE_ZSTD_SEARCH_DIR)
    file(GLOB _ZSTD_LIB "${_OPENEYE_ZSTD_SEARCH_DIR}/libzstd*.dylib" "${_OPENEYE_ZSTD_SEARCH_DIR}/libzstd*.so")
    if(_ZSTD_LIB)
        list(GET _ZSTD_LIB 0 OEZSTD_LIBRARY)
        message(STATUS "OpenEye: Found zstd: ${OEZSTD_LIBRARY}")
    endif()
endif()
if(NOT OEZSTD_LIBRARY)
    find_library(OEZSTD_LIBRARY
        NAMES zstd_static zstd
        PATHS
            ${_OPENEYE_LIB_SEARCH_PATHS}
            ${OPENEYE_ROOT}/lib
            ${OE_DIR}/lib
            $ENV{OPENEYE_ROOT}/lib
            $ENV{OE_DIR}/lib
            ${_OPENEYE_SYSTEM_LIB_FALLBACKS}
        NO_DEFAULT_PATH
    )
endif()

# Restore CMAKE_FIND_LIBRARY_SUFFIXES before finding system libraries
set(CMAKE_FIND_LIBRARY_SUFFIXES ${_SAVED_CMAKE_FIND_LIBRARY_SUFFIXES})

# Find system zlib. On Windows zlib isn't a system library, so fall back to
# FetchContent so downstream projects don't need to provide it themselves.
# Skip ZLIB in script mode (for tests) since find_package cannot create targets.
#
# Note: The OpenEye Windows SDK ships z.lib (a static archive of OpenEye's
# bundled zlib) but no zlib headers. Consumers like oemaestro include
# <zlib.h> directly, so we still fetch zlib v1.3.1 sources for headers and
# use the resulting zlibstatic target for linking. This keeps a single
# zlib copy in the final .pyd — matching what OE's own DLLs statically link.
if(NOT CMAKE_SCRIPT_MODE_FILE)
    find_package(ZLIB QUIET)
    if(NOT ZLIB_FOUND)
        if(WIN32)
            message(STATUS "OpenEye: ZLIB not found; fetching zlib v1.3.1 for Windows build")
            include(FetchContent)
            set(ZLIB_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
            set(SKIP_INSTALL_ALL ON CACHE BOOL "" FORCE)
            FetchContent_Declare(
                zlib
                GIT_REPOSITORY https://github.com/madler/zlib.git
                GIT_TAG v1.3.1
                GIT_SHALLOW TRUE
            )
            FetchContent_MakeAvailable(zlib)
            set(ZLIB_INCLUDE_DIR "${zlib_SOURCE_DIR};${zlib_BINARY_DIR}" CACHE PATH "" FORCE)
            set(ZLIB_LIBRARY zlibstatic CACHE STRING "" FORCE)
            set(ZLIB_FOUND TRUE CACHE BOOL "" FORCE)
            if(NOT TARGET ZLIB::ZLIB)
                add_library(ZLIB::ZLIB ALIAS zlibstatic)
            endif()
        else()
            find_package(ZLIB REQUIRED)
        endif()
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenEye
    REQUIRED_VARS
        OPENEYE_INCLUDE_DIR
        OECHEM_LIBRARY
        OESYSTEM_LIBRARY
        OEPLATFORM_LIBRARY
        OEMATH_LIBRARY
)

# Determine library type based on file extension. Windows .lib is always a
# static archive here (the SDK ships no DLL import libraries and the wheel
# ships no .lib files at all).
if(OpenEye_FOUND)
    get_filename_component(OECHEM_NAME "${OECHEM_LIBRARY}" NAME)
    if(OECHEM_NAME MATCHES "\\.dylib$" OR OECHEM_NAME MATCHES "\\.so$" OR OECHEM_NAME MATCHES "\\.so\\.")
        set(OPENEYE_LIBRARY_TYPE SHARED)
        message(STATUS "OpenEye: Using shared libraries (dynamic linking)")
    else()
        if(OPENEYE_USE_SHARED)
            message(FATAL_ERROR
                "OpenEye: OPENEYE_USE_SHARED is ON, but shared OpenEye libraries were not found. "
                "Resolved OEChem to ${OECHEM_LIBRARY}. Set OPENEYE_LIB_DIR or OPENEYE_RUNTIME_LIB_DIR "
                "to the openeye-toolkits shared library directory.")
        endif()
        set(OPENEYE_LIBRARY_TYPE STATIC)
        message(STATUS "OpenEye: Using static libraries")
    endif()

    # Extract OEChem version from library name (e.g., liboechem-4.3.0.1.dylib)
    # This matches what OEChemGetVersion() returns at runtime
    get_filename_component(OECHEM_NAME "${OECHEM_LIBRARY}" NAME)
    string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+" OpenEye_VERSION "${OECHEM_NAME}")
    if(NOT OpenEye_VERSION)
        # Try shorter version format (e.g., 4.3.0)
        string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" OpenEye_VERSION "${OECHEM_NAME}")
    endif()
    if(OpenEye_VERSION)
        message(STATUS "OpenEye: OEChem version ${OpenEye_VERSION}")
    else()
        # Try to extract from path as fallback
        string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" OpenEye_VERSION "${OPENEYE_INCLUDE_DIR}")
        if(OpenEye_VERSION)
            message(STATUS "OpenEye: Toolkit version ${OpenEye_VERSION} (from path)")
        endif()
    endif()

    # Export the library type for use in other CMake files, including
    # script-mode tests where imported targets are intentionally not created.
    set(OpenEye_LIBRARY_TYPE ${OPENEYE_LIBRARY_TYPE} CACHE STRING "OpenEye library type (SHARED or STATIC)")
endif()

# Detect the OpenEye SDK major year (e.g., 2024, 2025) and full release. The
# dep graph changes across SDK majors (notably OESpruce in 2025.2+), so
# downstream logic may condition on OpenEye_SDK_VERSION / OpenEye_SDK_MAJOR.
#
# Sources, in priority order:
#   1. README.txt at the SDK root ("OpenEye Toolkits vYYYY.N.N ..."). Shipped
#      in every SDK tarball and carries the full release triplet.
#   2. The install-path regex, accepting .../toolkits/YYYY.N.N/... through
#      .../toolkits/YYYY/... Covers Python-package-style install layouts.
#   3. A bare-year conservative fallback. Deliberately NOT ".N" so
#      VERSION_GREATER_EQUAL "2025.2" evaluates FALSE by default — callers who
#      need the newer dep graph must set the version explicitly rather than
#      silently inheriting it from a failed detection.
set(OpenEye_SDK_VERSION "")
set(OpenEye_SDK_MAJOR "")
if(OPENEYE_INCLUDE_DIR)
    # PATH_SUFFIXES openeye means find_path can land on either include/ (real
    # SDK) or include/openeye/ (synthetic test fixtures). Walk upward until we
    # find a directory containing README.txt, checking at most three levels.
    set(_OE_SDK_SEARCH "${OPENEYE_INCLUDE_DIR}")
    set(_OE_README "")
    foreach(_ _ _ _)
        get_filename_component(_OE_SDK_SEARCH "${_OE_SDK_SEARCH}/.." ABSOLUTE)
        if(EXISTS "${_OE_SDK_SEARCH}/README.txt")
            set(_OE_README "${_OE_SDK_SEARCH}/README.txt")
            break()
        endif()
        if(_OE_SDK_SEARCH STREQUAL "/")
            break()
        endif()
    endforeach()
    if(_OE_README)
        file(STRINGS "${_OE_README}" _OE_README_LINE
            REGEX "OpenEye Toolkits[ \t]+v[0-9]+\\.[0-9]+(\\.[0-9]+)?")
        if(_OE_README_LINE)
            string(REGEX MATCH "v([0-9]+\\.[0-9]+(\\.[0-9]+)?)" _MATCH "${_OE_README_LINE}")
            set(OpenEye_SDK_VERSION "${CMAKE_MATCH_1}")
            string(REGEX MATCH "^([0-9]+)" _MAJOR_MATCH "${OpenEye_SDK_VERSION}")
            set(OpenEye_SDK_MAJOR "${CMAKE_MATCH_1}")
        endif()
    endif()
endif()
# Install-path fallback. Accepts .../toolkits/2025.2.1/..., .../toolkits/2025.2/...,
# or bare-year .../toolkits/2025/... The minor/patch components are optional so a
# year-only install path still populates OpenEye_SDK_MAJOR (the dep-graph gate
# only needs the major year).
if(NOT OpenEye_SDK_MAJOR AND OPENEYE_INCLUDE_DIR)
    string(REGEX MATCH "/(20[0-9][0-9](\\.[0-9]+(\\.[0-9]+)?)?)" _MATCH "${OPENEYE_INCLUDE_DIR}")
    if(CMAKE_MATCH_1)
        set(OpenEye_SDK_VERSION "${CMAKE_MATCH_1}")
        string(REGEX MATCH "^([0-9]+)" _MAJOR_MATCH "${OpenEye_SDK_VERSION}")
        set(OpenEye_SDK_MAJOR "${CMAKE_MATCH_1}")
    endif()
endif()
if(NOT OpenEye_SDK_MAJOR)
    set(OpenEye_SDK_MAJOR "2025")
    set(OpenEye_SDK_VERSION "2025")
    message(WARNING "OpenEye: Could not detect SDK version from ${OPENEYE_INCLUDE_DIR}; "
        "defaulting to OpenEye_SDK_MAJOR=${OpenEye_SDK_MAJOR} and "
        "OpenEye_SDK_VERSION=${OpenEye_SDK_VERSION} (falls back to the pre-2025.2 dep graph)")
else()
    message(STATUS "OpenEye: Detected SDK version ${OpenEye_SDK_VERSION}")
endif()

if(OpenEye_FOUND AND NOT TARGET OpenEye::OEChem AND NOT CMAKE_SCRIPT_MODE_FILE)
    # Create imported target for zstd if found. The bundled static zstd uses
    # pthread_create/join internally (see libzstd pool.c), so pre-glibc-2.34
    # systems (RHEL 8, Ubuntu 20.04) need an explicit -lpthread on the link
    # line. Pull in Threads::Threads so consumers of OpenEye::zstd do not need
    # to know about this.
    if(OEZSTD_LIBRARY AND NOT TARGET OpenEye::zstd)
        find_package(Threads)
        add_library(OpenEye::zstd UNKNOWN IMPORTED)
        set_target_properties(OpenEye::zstd PROPERTIES
            IMPORTED_LOCATION "${OEZSTD_LIBRARY}"
        )
        if(TARGET Threads::Threads)
            set_property(TARGET OpenEye::zstd APPEND PROPERTY
                INTERFACE_LINK_LIBRARIES Threads::Threads
            )
        endif()
    endif()

    # OEPlatform depends on zlib and zstd. Even on Windows (where oeplatform.lib
    # is a static archive with unresolved gzstd_*/zlib externs), we must link
    # both. The OpenEye Windows SDK ships z.lib and zstd_static.lib precisely
    # so downstream consumers can satisfy those externs.
    add_library(OpenEye::OEPlatform UNKNOWN IMPORTED)
    set_target_properties(OpenEye::OEPlatform PROPERTIES
        IMPORTED_LOCATION "${OEPLATFORM_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
    )
    if(OEZSTD_LIBRARY)
        set_property(TARGET OpenEye::OEPlatform APPEND PROPERTY
            INTERFACE_LINK_LIBRARIES "OpenEye::zstd;ZLIB::ZLIB"
        )
    else()
        set_property(TARGET OpenEye::OEPlatform APPEND PROPERTY
            INTERFACE_LINK_LIBRARIES "ZLIB::ZLIB"
        )
    endif()

    # OEPlatform's Windows hostinfo uses Winsock and Netbios.
    if(WIN32)
        set_property(TARGET OpenEye::OEPlatform APPEND PROPERTY
            INTERFACE_LINK_LIBRARIES "ws2_32;netapi32"
        )
    endif()

    add_library(OpenEye::OESystem UNKNOWN IMPORTED)
    set_target_properties(OpenEye::OESystem PROPERTIES
        IMPORTED_LOCATION "${OESYSTEM_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "OpenEye::OEPlatform"
    )

    # OEMath depends on OESystem (oemath/matrix.h includes oesystem.h)
    add_library(OpenEye::OEMath UNKNOWN IMPORTED)
    set_target_properties(OpenEye::OEMath PROPERTIES
        IMPORTED_LOCATION "${OEMATH_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "OpenEye::OESystem"
    )

    # OEChem depends on OESystem and OEMath
    add_library(OpenEye::OEChem UNKNOWN IMPORTED)
    set_target_properties(OpenEye::OEChem PROPERTIES
        IMPORTED_LOCATION "${OECHEM_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "OpenEye::OESystem;OpenEye::OEMath"
    )

    # Optional: OEGraphSim depends on OEChem
    if(OEGRAPHSIM_LIBRARY)
        add_library(OpenEye::OEGraphSim UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEGraphSim PROPERTIES
            IMPORTED_LOCATION "${OEGRAPHSIM_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem"
        )
        set(OpenEye_GraphSim_FOUND TRUE)
    endif()

    # Optional: OEMedChem depends on OEChem
    if(OEMEDCHEM_LIBRARY)
        add_library(OpenEye::OEMedChem UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEMedChem PROPERTIES
            IMPORTED_LOCATION "${OEMEDCHEM_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem"
        )
        set(OpenEye_MedChem_FOUND TRUE)
    endif()

    # Optional: OEFizzChem depends on OEChem
    if(OEFIZZCHEM_LIBRARY)
        add_library(OpenEye::OEFizzChem UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEFizzChem PROPERTIES
            IMPORTED_LOCATION "${OEFIZZCHEM_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem"
        )
    endif()

    # Optional: OEGrid depends on OESystem and OEFizzChem (if found)
    if(OEGRID_LIBRARY)
        add_library(OpenEye::OEGrid UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEGrid PROPERTIES
            IMPORTED_LOCATION "${OEGRID_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
        )
        if(OEFIZZCHEM_LIBRARY)
            set_property(TARGET OpenEye::OEGrid PROPERTY
                INTERFACE_LINK_LIBRARIES "OpenEye::OESystem;OpenEye::OEFizzChem"
            )
        else()
            set_property(TARGET OpenEye::OEGrid PROPERTY
                INTERFACE_LINK_LIBRARIES "OpenEye::OESystem"
            )
        endif()
        set(OpenEye_Grid_FOUND TRUE)
    endif()

    # Optional: OEBio depends on OEChem and OEGrid (uses OESkewGrid, OEScalarGrid, OEXtal)
    if(OEBIO_LIBRARY)
        add_library(OpenEye::OEBio UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEBio PROPERTIES
            IMPORTED_LOCATION "${OEBIO_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
        )
        if(OEGRID_LIBRARY)
            set_property(TARGET OpenEye::OEBio PROPERTY
                INTERFACE_LINK_LIBRARIES "OpenEye::OEChem;OpenEye::OEGrid"
            )
        else()
            set_property(TARGET OpenEye::OEBio PROPERTY
                INTERFACE_LINK_LIBRARIES "OpenEye::OEChem"
            )
        endif()
        set(OpenEye_Bio_FOUND TRUE)
    endif()

    # v1.1.0: Geometry and optimization targets
    if(OEOPT_LIBRARY)
        add_library(OpenEye::OEOpt UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEOpt PROPERTIES
            IMPORTED_LOCATION "${OEOPT_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OESystem"
        )
        set(OpenEye_Opt_FOUND TRUE)
    endif()

    if(OEMOLPOTENTIAL_LIBRARY AND TARGET OpenEye::OEOpt)
        add_library(OpenEye::OEMolPotential UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEMolPotential PROPERTIES
            IMPORTED_LOCATION "${OEMOLPOTENTIAL_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem;OpenEye::OEOpt"
        )
        set(OpenEye_MolPotential_FOUND TRUE)
    endif()

    if(OEHERMITE_LIBRARY AND TARGET OpenEye::OEOpt)
        add_library(OpenEye::OEHermite UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEHermite PROPERTIES
            IMPORTED_LOCATION "${OEHERMITE_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem;OpenEye::OEOpt"
        )
        set(OpenEye_Hermite_FOUND TRUE)
    endif()

    # OEShape's umbrella header includes oebio.h (via sitehopperdatabase_base.h),
    # which instantiates OEFieldType<OEBio::OEDesignUnit>. Users of OEShape must
    # therefore link OEBio for the vtable symbol, so OEBio is a hard dep here.
    if(OESHAPE_LIBRARY
            AND TARGET OpenEye::OEBio
            AND TARGET OpenEye::OEGrid
            AND TARGET OpenEye::OEOpt
            AND TARGET OpenEye::OEMolPotential
            AND TARGET OpenEye::OEHermite)
        add_library(OpenEye::OEShape UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEShape PROPERTIES
            IMPORTED_LOCATION "${OESHAPE_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem;OpenEye::OEBio;OpenEye::OEGrid;OpenEye::OEOpt;OpenEye::OEMolPotential;OpenEye::OEHermite"
        )
        set(OpenEye_Shape_FOUND TRUE)
    endif()

    if(OEZAP_LIBRARY AND TARGET OpenEye::OEGrid)
        add_library(OpenEye::OEZap UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEZap PROPERTIES
            IMPORTED_LOCATION "${OEZAP_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem;OpenEye::OEGrid"
        )
        set(OpenEye_Zap_FOUND TRUE)
    endif()

    # OESpicoli references OEBio symbols (OEGetResidues, OEIsNTerminalAtom, OEIsWater)
    # in its archive, so OEBio is a hard runtime dep — guard on it here.
    if(OESPICOLI_LIBRARY AND TARGET OpenEye::OEZap AND TARGET OpenEye::OEBio)
        add_library(OpenEye::OESpicoli UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OESpicoli PROPERTIES
            IMPORTED_LOCATION "${OESPICOLI_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem;OpenEye::OEZap;OpenEye::OEBio"
        )
        set(OpenEye_Spicoli_FOUND TRUE)
    endif()

    if(OESITEHOPPER_LIBRARY AND TARGET OpenEye::OEShape AND TARGET OpenEye::OESpicoli)
        add_library(OpenEye::OESiteHopper UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OESiteHopper PROPERTIES
            IMPORTED_LOCATION "${OESITEHOPPER_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem;OpenEye::OEShape;OpenEye::OESpicoli;${CMAKE_DL_LIBS}"
        )
        set(OpenEye_SiteHopper_FOUND TRUE)
    endif()

    # v1.1.0: Force-field / conformer / protonation targets
    # liboemmff.a references 175 undefined OEMolPotential::* symbols
    # (e.g. OEForceField::OEForceField(), OEGenericFF2::AddMolFunc), so
    # OEMolPotential is a hard runtime dep and appears in INTERFACE_LINK_LIBRARIES.
    # Not guarded here because OEMolPotential is a core lib that should always
    # be present alongside OEMMFF in any SDK shipping force-field support.
    if(OEMMFF_LIBRARY AND TARGET OpenEye::OEMolPotential)
        add_library(OpenEye::OEMMFF UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEMMFF PROPERTIES
            IMPORTED_LOCATION "${OEMMFF_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem;OpenEye::OEMolPotential"
        )
        set(OpenEye_MMFF_FOUND TRUE)
    endif()

    # Newer force-field stacks split SMIRNOFF, Amber, and AM1-BCC support into
    # separate archives. Model them explicitly so static consumers do not need
    # to know private toolkit-library closure details.
    if(OEAMBER_LIBRARY AND TARGET OpenEye::OEMolPotential)
        add_library(OpenEye::OEAmber UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEAmber PROPERTIES
            IMPORTED_LOCATION "${OEAMBER_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem;OpenEye::OEMolPotential"
        )
        set(OpenEye_Amber_FOUND TRUE)
    endif()

    if(OESMIRNOFF_LIBRARY AND TARGET OpenEye::OEMolPotential AND TARGET OpenEye::OEAmber)
        add_library(OpenEye::OESmirnoff UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OESmirnoff PROPERTIES
            IMPORTED_LOCATION "${OESMIRNOFF_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem;OpenEye::OEMolPotential;OpenEye::OEAmber"
        )
        set(OpenEye_Smirnoff_FOUND TRUE)
    endif()

    if(OEAM1_LIBRARY AND TARGET OpenEye::OEMolPotential)
        add_library(OpenEye::OEAM1 UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEAM1 PROPERTIES
            IMPORTED_LOCATION "${OEAM1_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem;OpenEye::OEMolPotential"
        )
        set(OpenEye_AM1_FOUND TRUE)
    endif()

    if(OEAM1BCC_LIBRARY AND TARGET OpenEye::OEAM1)
        add_library(OpenEye::OEAM1BCC UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEAM1BCC PROPERTIES
            IMPORTED_LOCATION "${OEAM1BCC_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem;OpenEye::OEAM1"
        )
        set(OpenEye_AM1BCC_FOUND TRUE)
    endif()

    if(OEFF_LIBRARY AND TARGET OpenEye::OEMolPotential AND TARGET OpenEye::OESmirnoff AND TARGET OpenEye::OEAmber)
        add_library(OpenEye::OEFF UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEFF PROPERTIES
            IMPORTED_LOCATION "${OEFF_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem;OpenEye::OEMolPotential;OpenEye::OESmirnoff;OpenEye::OEAmber"
        )
        set(OpenEye_FF_FOUND TRUE)
    endif()

    # liboeszybki.a references 11 undefined OEBio symbols
    # (e.g. OEBio::OEDesignUnitImpl::SetComponentsFromData, OEBio::OEAtomMatchResidue)
    # in its archive, so OEBio is a hard runtime dep — guard on it here.
    if(OESZYBKI_LIBRARY AND TARGET OpenEye::OEMMFF AND TARGET OpenEye::OEFF AND TARGET OpenEye::OEBio)
        add_library(OpenEye::OESzybki UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OESzybki PROPERTIES
            IMPORTED_LOCATION "${OESZYBKI_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem;OpenEye::OEMMFF;OpenEye::OEFF;OpenEye::OEBio"
        )
        set(OpenEye_Szybki_FOUND TRUE)
    endif()

    if(OEQUACPAC_LIBRARY AND TARGET OpenEye::OESzybki AND TARGET OpenEye::OEAmber AND TARGET OpenEye::OEAM1BCC)
        add_library(OpenEye::OEQuacpac UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEQuacpac PROPERTIES
            IMPORTED_LOCATION "${OEQUACPAC_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem;OpenEye::OESzybki;OpenEye::OEAmber;OpenEye::OEAM1BCC"
        )
        set(OpenEye_Quacpac_FOUND TRUE)
    endif()

    if(OEOMEGA2_LIBRARY AND TARGET OpenEye::OEMMFF AND TARGET OpenEye::OESmirnoff)
        add_library(OpenEye::OEOmega2 UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEOmega2 PROPERTIES
            IMPORTED_LOCATION "${OEOMEGA2_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem;OpenEye::OEMMFF;OpenEye::OESmirnoff"
        )
        set(OpenEye_Omega2_FOUND TRUE)
    endif()

    # liboesheffield.a references undefined symbols from four separate libs:
    #   * OEFizzChem provides OEFizzChem::OEDerefGridHandle, OEReleaseGridHandle
    #   * OEGrid provides OESystem::oe_read_grid_object and
    #     oe_convert_grid_object_to_grid_structure. Crucially these are NOT
    #     transitive through OEFizzChem — liboefizzchem.a has zero OESystem::oe_*
    #     grid undefs, so OEGrid must be a direct dep. Do not drop it.
    #   * OEZap provides OEPB::oe_make_zap, oe_make_area
    #   * OEMolPotential provides 25 undefined OEMolPotential::* refs
    #     (e.g. OEMolPotential::OEMolFunc::SetVerbose)
    # All four are hard runtime deps; OEFizzChem/OEGrid/OEZap are guarded above
    # and OEMolPotential is included in INTERFACE_LINK_LIBRARIES.
    if(OESHEFFIELD_LIBRARY
            AND TARGET OpenEye::OEMolPotential
            AND TARGET OpenEye::OEFizzChem
            AND TARGET OpenEye::OEGrid
            AND TARGET OpenEye::OEZap)
        add_library(OpenEye::OESheffield UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OESheffield PROPERTIES
            IMPORTED_LOCATION "${OESHEFFIELD_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem;OpenEye::OEMolPotential;OpenEye::OEFizzChem;OpenEye::OEGrid;OpenEye::OEZap"
        )
        set(OpenEye_Sheffield_FOUND TRUE)
    endif()

    # OESpruce's transitive deps expanded in SDK 2025.2. Older SDKs only need
    # OEChem + OEBio + OESiteHopper; 2025.2+ also pulls OEQuacpac, OEMMFF,
    # OEOmega2, OESheffield via builder/designunit code paths.
    if(OESPRUCE_LIBRARY)
        set(_OPENEYE_SPRUCE_DEPS_AVAILABLE FALSE)
        if(OpenEye_SDK_VERSION VERSION_GREATER_EQUAL "2025.2")
            if(TARGET OpenEye::OEBio
                    AND TARGET OpenEye::OESiteHopper
                    AND TARGET OpenEye::OEQuacpac
                    AND TARGET OpenEye::OEMMFF
                    AND TARGET OpenEye::OEOmega2
                    AND TARGET OpenEye::OESheffield)
                set(_OPENEYE_SPRUCE_DEPS_AVAILABLE TRUE)
            endif()
        else()
            if(TARGET OpenEye::OEBio AND TARGET OpenEye::OESiteHopper)
                set(_OPENEYE_SPRUCE_DEPS_AVAILABLE TRUE)
            endif()
        endif()
    endif()

    if(OESPRUCE_LIBRARY AND _OPENEYE_SPRUCE_DEPS_AVAILABLE)
        add_library(OpenEye::OESpruce UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OESpruce PROPERTIES
            IMPORTED_LOCATION "${OESPRUCE_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
        )
        if(OpenEye_SDK_VERSION VERSION_GREATER_EQUAL "2025.2")
            set_property(TARGET OpenEye::OESpruce PROPERTY INTERFACE_LINK_LIBRARIES
                "OpenEye::OEChem;OpenEye::OEBio;OpenEye::OESiteHopper;OpenEye::OEQuacpac;OpenEye::OEMMFF;OpenEye::OEOmega2;OpenEye::OESheffield"
            )
        else()
            set_property(TARGET OpenEye::OESpruce PROPERTY INTERFACE_LINK_LIBRARIES
                "OpenEye::OEChem;OpenEye::OEBio;OpenEye::OESiteHopper"
            )
        endif()
        set(OpenEye_Spruce_FOUND TRUE)
    endif()

    # v1.1.0: Depiction / nomenclature targets
    if(OEDEPICT_LIBRARY)
        add_library(OpenEye::OEDepict UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEDepict PROPERTIES
            IMPORTED_LOCATION "${OEDEPICT_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem"
        )
        set(OpenEye_Depict_FOUND TRUE)
    endif()

    if(OEIUPAC_LIBRARY)
        add_library(OpenEye::OEIUPAC UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEIUPAC PROPERTIES
            IMPORTED_LOCATION "${OEIUPAC_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem"
        )
        set(OpenEye_IUPAC_FOUND TRUE)
    endif()

endif()

mark_as_advanced(
    OPENEYE_INCLUDE_DIR
    OECHEM_LIBRARY
    OESYSTEM_LIBRARY
    OEPLATFORM_LIBRARY
    OEGRAPHSIM_LIBRARY
    OEMEDCHEM_LIBRARY
    OEBIO_LIBRARY
    OEGRID_LIBRARY
    OEFIZZCHEM_LIBRARY
    OEMATH_LIBRARY
    OEZSTD_LIBRARY
    OEOPT_LIBRARY
    OEMOLPOTENTIAL_LIBRARY
    OEHERMITE_LIBRARY
    OESHAPE_LIBRARY
    OEZAP_LIBRARY
    OESPICOLI_LIBRARY
    OESITEHOPPER_LIBRARY
    OEMMFF_LIBRARY
    OEFF_LIBRARY
    OESZYBKI_LIBRARY
    OEQUACPAC_LIBRARY
    OEOMEGA2_LIBRARY
    OESHEFFIELD_LIBRARY
    OESPRUCE_LIBRARY
    OESMIRNOFF_LIBRARY
    OEAMBER_LIBRARY
    OEAM1_LIBRARY
    OEAM1BCC_LIBRARY
    OEDEPICT_LIBRARY
    OEIUPAC_LIBRARY
)
