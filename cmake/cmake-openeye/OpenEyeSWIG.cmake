# OpenEyeSWIG.cmake
# Provides openeye_add_swig_module() for building SWIG Python bindings
# against OpenEye C++ toolkits.
#
# This module encapsulates all SWIG boilerplate: finding SWIG/Python,
# stable ABI support, platform-specific linking, RPATH configuration,
# post-build copies for editable installs, install rules, and
# _build_info.py generation.

#[============================================================================[
  openeye_add_swig_module(
      NAME <name>
      SWIG_FILE <file>
      LINK_LIBS <libs...>
      PYTHON_OUTPUT_DIR <dir>
      [STABLE_ABI <bool>]
      [SWIG_FLAGS <flags...>]
      [COMPILE_DEFS <defs...>]
      [EXTRA_INSTALL_TARGETS <targets...>]
      [EXPECTED_LIB_VARS <vars...>]
      [INIT_PY <path>]
  )
#]============================================================================]

function(openeye_add_swig_module)
    # Parse arguments
    set(_OPTIONS "")
    set(_ONE_VALUE_ARGS NAME SWIG_FILE PYTHON_OUTPUT_DIR STABLE_ABI INIT_PY)
    set(_MULTI_VALUE_ARGS LINK_LIBS SWIG_FLAGS COMPILE_DEFS EXTRA_INSTALL_TARGETS EXPECTED_LIB_VARS)
    cmake_parse_arguments(ARG "${_OPTIONS}" "${_ONE_VALUE_ARGS}" "${_MULTI_VALUE_ARGS}" ${ARGN})

    # Validate required arguments
    if(NOT ARG_NAME)
        message(FATAL_ERROR "openeye_add_swig_module: NAME is required")
    endif()
    if(NOT ARG_SWIG_FILE)
        message(FATAL_ERROR "openeye_add_swig_module: SWIG_FILE is required")
    endif()
    if(NOT ARG_LINK_LIBS)
        message(FATAL_ERROR "openeye_add_swig_module: LINK_LIBS is required")
    endif()
    if(NOT ARG_PYTHON_OUTPUT_DIR)
        message(FATAL_ERROR "openeye_add_swig_module: PYTHON_OUTPUT_DIR is required")
    endif()

    # Default STABLE_ABI to OFF if not provided
    if(NOT DEFINED ARG_STABLE_ABI)
        set(ARG_STABLE_ABI OFF)
    endif()

    set(_TARGET_NAME ${ARG_NAME}_python)
    set(_MODULE_NAME _${ARG_NAME})

    # -------------------------------------------------------------------------
    # 1. Find SWIG
    # -------------------------------------------------------------------------
    if(ARG_STABLE_ABI)
        find_package(SWIG 4.2 REQUIRED)
    else()
        find_package(SWIG 4.0 REQUIRED)
    endif()

    # -------------------------------------------------------------------------
    # 2. Find Python3 (with SABIModule when STABLE_ABI is ON)
    # -------------------------------------------------------------------------
    # STABLE_ABI requires Python3::SABIModule, introduced in CMake 3.26.
    # On Windows this is load-bearing: without it, the extension links
    # python3XX.lib (versioned) instead of python3.lib (stable-ABI), which
    # silently breaks the abi3 wheel contract.
    if(ARG_STABLE_ABI AND CMAKE_VERSION VERSION_LESS "3.26")
        message(FATAL_ERROR
            "OpenEyeSWIG: STABLE_ABI=ON requires CMake 3.26+ for the "
            "Python3::SABIModule imported target. CMake in use: "
            "${CMAKE_VERSION}. Either upgrade CMake or set STABLE_ABI=OFF.")
    endif()
    if(ARG_STABLE_ABI)
        find_package(Python3 COMPONENTS Interpreter Development.SABIModule REQUIRED)
    else()
        find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
    endif()

    include(${SWIG_USE_FILE})

    # -------------------------------------------------------------------------
    # 3. SWIG source properties
    # -------------------------------------------------------------------------
    set_source_files_properties(${ARG_SWIG_FILE} PROPERTIES
        CPLUSPLUS ON
        SWIG_MODULE_NAME ${_MODULE_NAME}
    )

    # Set SWIG flags (caller provides all flags via SWIG_FLAGS)
    set(CMAKE_SWIG_FLAGS "")
    if(ARG_SWIG_FLAGS)
        set(CMAKE_SWIG_FLAGS ${ARG_SWIG_FLAGS})
    endif()

    # -------------------------------------------------------------------------
    # 4. Create SWIG library (MODULE type)
    # -------------------------------------------------------------------------
    swig_add_library(${_TARGET_NAME}
        TYPE MODULE
        LANGUAGE python
        OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}
        SOURCES ${ARG_SWIG_FILE}
    )

    # -------------------------------------------------------------------------
    # 5. Stable ABI (Py_LIMITED_API)
    # -------------------------------------------------------------------------
    if(ARG_STABLE_ABI)
        target_compile_definitions(${_TARGET_NAME} PRIVATE
            Py_LIMITED_API=0x030A0000
        )
        message(STATUS "Building ${ARG_NAME} with Python stable ABI (abi3) for Python 3.10+")
    endif()

    # -------------------------------------------------------------------------
    # 6. Include directories
    # -------------------------------------------------------------------------
    target_include_directories(${_TARGET_NAME} PRIVATE
        ${Python3_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/include
    )

    # -------------------------------------------------------------------------
    # 7. Extra compile definitions
    # -------------------------------------------------------------------------
    if(ARG_COMPILE_DEFS)
        target_compile_definitions(${_TARGET_NAME} PRIVATE ${ARG_COMPILE_DEFS})
    endif()

    # -------------------------------------------------------------------------
    # 8. Platform-specific linking
    # -------------------------------------------------------------------------
    if(APPLE)
        # On macOS, use -undefined dynamic_lookup to allow Python symbols to be
        # resolved at runtime by the Python interpreter that loads the module.
        # This is critical for Python 3.13+ multi-phase initialization.
        target_link_libraries(${_TARGET_NAME}
            PRIVATE
                ${ARG_LINK_LIBS}
        )
        target_link_options(${_TARGET_NAME} PRIVATE
            -undefined dynamic_lookup
        )
    elseif(ARG_STABLE_ABI)
        # Py_LIMITED_API is defined on the target; link Python3::SABIModule
        # so the extension resolves against python3.lib -> python3.dll. Linking
        # Python3::Python here would drag in python3XX.lib (versioned) and
        # silently break the abi3 wheel contract on Windows.
        target_link_libraries(${_TARGET_NAME}
            PRIVATE
                ${ARG_LINK_LIBS}
                Python3::SABIModule
        )
    else()
        target_link_libraries(${_TARGET_NAME}
            PRIVATE
                ${ARG_LINK_LIBS}
                Python3::Python
        )
    endif()

    # -------------------------------------------------------------------------
    # 9. Output properties: prefix, suffix, output name
    # -------------------------------------------------------------------------
    set_target_properties(${_TARGET_NAME} PROPERTIES
        OUTPUT_NAME "${_MODULE_NAME}"
        PREFIX ""
    )
    if(APPLE)
        set_target_properties(${_TARGET_NAME} PROPERTIES SUFFIX ".so")
    elseif(WIN32)
        set_target_properties(${_TARGET_NAME} PROPERTIES SUFFIX ".pyd")
    endif()

    # -------------------------------------------------------------------------
    # 10. RPATH configuration
    # -------------------------------------------------------------------------
    # Use OpenEyePython_PLATFORM if available, otherwise query Python at configure time
    if(NOT DEFINED _OE_SWIG_PLATFORM)
        if(DEFINED OpenEyePython_PLATFORM AND OpenEyePython_PLATFORM)
            set(_OE_SWIG_PLATFORM "${OpenEyePython_PLATFORM}")
        elseif(DEFINED OPENEYE_RUNTIME_LIB_DIR AND OPENEYE_RUNTIME_LIB_DIR)
            get_filename_component(_OE_SWIG_PLATFORM "${OPENEYE_RUNTIME_LIB_DIR}" NAME)
        elseif(OpenEye_LIBRARY_TYPE STREQUAL "SHARED")
            execute_process(
                COMMAND ${Python3_EXECUTABLE} -c
                    "from openeye import libs; import os; print(os.path.basename(libs.FindOpenEyeDLLSDirectory()))"
                OUTPUT_VARIABLE _OE_SWIG_PLATFORM
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
            )
        endif()
    endif()

    if(APPLE)
        if(OpenEye_LIBRARY_TYPE STREQUAL "SHARED")
            set_target_properties(${_TARGET_NAME} PROPERTIES
                INSTALL_RPATH "@loader_path;@loader_path/../openeye/libs/${_OE_SWIG_PLATFORM}"
                BUILD_WITH_INSTALL_RPATH TRUE
            )
            message(STATUS "${ARG_NAME}: RPATH set for OpenEye Python package: @loader_path/../openeye/libs/${_OE_SWIG_PLATFORM}")
        else()
            set_target_properties(${_TARGET_NAME} PROPERTIES
                INSTALL_RPATH "@loader_path"
                BUILD_WITH_INSTALL_RPATH TRUE
            )
        endif()
    elseif(WIN32)
        # Windows has no RPATH. DLL resolution happens at Python import time via
        # openeye.libs' os.add_dll_directory() side effect — the consumer's
        # python/<pkg>/__init__.py must `import openeye.libs` before loading the
        # SWIG .pyd module.
        message(STATUS "${ARG_NAME}: Windows — no RPATH; DLL resolution via openeye.libs at import time")
    elseif(UNIX)
        if(OpenEye_LIBRARY_TYPE STREQUAL "SHARED")
            set_target_properties(${_TARGET_NAME} PROPERTIES
                INSTALL_RPATH "$ORIGIN:$ORIGIN/../openeye/libs/${_OE_SWIG_PLATFORM}"
                BUILD_WITH_INSTALL_RPATH TRUE
            )
        else()
            set_target_properties(${_TARGET_NAME} PROPERTIES
                INSTALL_RPATH "$ORIGIN"
                BUILD_WITH_INSTALL_RPATH TRUE
            )
        endif()
    endif()

    # -------------------------------------------------------------------------
    # 12. Post-build copy to PYTHON_OUTPUT_DIR
    # -------------------------------------------------------------------------
    add_custom_command(TARGET ${_TARGET_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory ${ARG_PYTHON_OUTPUT_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_BINARY_DIR}/${_MODULE_NAME}.py
            ${ARG_PYTHON_OUTPUT_DIR}/${ARG_NAME}.py
        COMMAND ${CMAKE_COMMAND} -E copy
            $<TARGET_FILE:${_TARGET_NAME}>
            ${ARG_PYTHON_OUTPUT_DIR}/$<TARGET_FILE_NAME:${_TARGET_NAME}>
        COMMENT "Copying ${ARG_NAME} Python module to ${ARG_PYTHON_OUTPUT_DIR}"
    )

    # Copy the static library for the main project lib (first LINK_LIB assumed to be the project lib)
    list(GET ARG_LINK_LIBS 0 _PRIMARY_LIB)
    if(TARGET ${_PRIMARY_LIB})
        add_custom_command(TARGET ${_TARGET_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy
                $<TARGET_FILE:${_PRIMARY_LIB}>
                ${ARG_PYTHON_OUTPUT_DIR}/$<TARGET_FILE_NAME:${_PRIMARY_LIB}>
            COMMENT "Copying lib${_PRIMARY_LIB} to ${ARG_PYTHON_OUTPUT_DIR}"
        )
    endif()

    # -------------------------------------------------------------------------
    # 13. EXTRA_INSTALL_TARGETS: copy alongside extension during post-build
    # -------------------------------------------------------------------------
    if(ARG_EXTRA_INSTALL_TARGETS)
        foreach(_EXTRA_TARGET ${ARG_EXTRA_INSTALL_TARGETS})
            if(TARGET ${_EXTRA_TARGET})
                add_custom_command(TARGET ${_TARGET_NAME} POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy
                        $<TARGET_FILE:${_EXTRA_TARGET}>
                        ${ARG_PYTHON_OUTPUT_DIR}/$<TARGET_FILE_NAME:${_EXTRA_TARGET}>
                    COMMENT "Copying ${_EXTRA_TARGET} to ${ARG_PYTHON_OUTPUT_DIR}"
                )
            endif()
        endforeach()
    endif()

    # -------------------------------------------------------------------------
    # 11. Editable dev RPATH fix (on PYTHON_OUTPUT_DIR copy only)
    # -------------------------------------------------------------------------
    # Determine absolute OpenEye lib dir for dev installs
    if(OpenEye_LIBRARY_TYPE STREQUAL "SHARED" AND _OE_SWIG_PLATFORM)
        if(DEFINED OPENEYE_RUNTIME_LIB_DIR AND OPENEYE_RUNTIME_LIB_DIR)
            set(_DEV_OE_LIB_DIR "${OPENEYE_RUNTIME_LIB_DIR}")
        elseif(DEFINED OpenEyePython_LIB_DIR AND OpenEyePython_LIB_DIR)
            set(_DEV_OE_LIB_DIR "${OpenEyePython_LIB_DIR}")
        else()
            execute_process(
                COMMAND ${Python3_EXECUTABLE} -c
                    "from openeye import libs; print(libs.FindOpenEyeDLLSDirectory())"
                OUTPUT_VARIABLE _DEV_OE_LIB_DIR
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
            )
        endif()

        if(_DEV_OE_LIB_DIR)
            if(APPLE)
                find_program(INSTALL_NAME_TOOL install_name_tool)
                if(INSTALL_NAME_TOOL)
                    add_custom_command(TARGET ${_TARGET_NAME} POST_BUILD
                        COMMAND ${INSTALL_NAME_TOOL} -add_rpath "${_DEV_OE_LIB_DIR}"
                            ${ARG_PYTHON_OUTPUT_DIR}/$<TARGET_FILE_NAME:${_TARGET_NAME}>
                        COMMAND codesign -f -s - ${ARG_PYTHON_OUTPUT_DIR}/$<TARGET_FILE_NAME:${_TARGET_NAME}>
                        COMMENT "Adding OpenEye RPATH for editable dev install (${ARG_NAME})"
                    )
                endif()
            elseif(UNIX)
                find_program(PATCHELF patchelf)
                if(PATCHELF)
                    add_custom_command(TARGET ${_TARGET_NAME} POST_BUILD
                        COMMAND ${PATCHELF} --add-rpath "${_DEV_OE_LIB_DIR}"
                            ${ARG_PYTHON_OUTPUT_DIR}/$<TARGET_FILE_NAME:${_TARGET_NAME}>
                        COMMENT "Adding OpenEye RPATH for editable dev install (${ARG_NAME})"
                    )
                endif()
            endif()
        endif()
    endif()

    # -------------------------------------------------------------------------
    # 14. Install rules
    # -------------------------------------------------------------------------
    # Install the SWIG extension module and primary library
    install(TARGETS ${_TARGET_NAME}
        LIBRARY DESTINATION ${ARG_NAME}
        RUNTIME DESTINATION ${ARG_NAME}
    )

    # Install the primary library target alongside the extension
    if(TARGET ${_PRIMARY_LIB})
        install(TARGETS ${_PRIMARY_LIB}
            LIBRARY DESTINATION ${ARG_NAME}
            RUNTIME DESTINATION ${ARG_NAME}
        )
    endif()

    # Install extra targets alongside the extension
    if(ARG_EXTRA_INSTALL_TARGETS)
        foreach(_EXTRA_TARGET ${ARG_EXTRA_INSTALL_TARGETS})
            if(TARGET ${_EXTRA_TARGET})
                install(TARGETS ${_EXTRA_TARGET}
                    LIBRARY DESTINATION ${ARG_NAME}
                    RUNTIME DESTINATION ${ARG_NAME}
                )
            endif()
        endforeach()
    endif()

    # Install the SWIG-generated Python file (renamed from _<name>.py to <name>.py)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${_MODULE_NAME}.py
        DESTINATION ${ARG_NAME}
        RENAME ${ARG_NAME}.py
    )

    # Install __init__.py if provided
    if(ARG_INIT_PY)
        install(FILES ${ARG_INIT_PY}
            DESTINATION ${ARG_NAME}
        )
    endif()

    # -------------------------------------------------------------------------
    # 15. _build_info.py generation
    # -------------------------------------------------------------------------
    # Determine version: prefer OPENEYE_TOOLKITS_VERSION (marketing), fall back to OpenEye_VERSION (library)
    if(DEFINED OPENEYE_TOOLKITS_VERSION AND OPENEYE_TOOLKITS_VERSION)
        set(_BUILD_INFO_VERSION "${OPENEYE_TOOLKITS_VERSION}")
    else()
        set(_BUILD_INFO_VERSION "${OpenEye_VERSION}")
    endif()

    # Collect expected library filenames for runtime compatibility (shared only).
    # Auto-detect all OpenEye library variables found by FindOpenEye.cmake,
    # then append any additional variables specified via EXPECTED_LIB_VARS.
    set(_ALL_OE_LIB_VARS
        OECHEM_LIBRARY OESYSTEM_LIBRARY OEPLATFORM_LIBRARY OEMATH_LIBRARY
        OEZSTD_LIBRARY OEGRAPHSIM_LIBRARY OEMEDCHEM_LIBRARY
        OEBIO_LIBRARY OEGRID_LIBRARY OEFIZZCHEM_LIBRARY
    )
    if(ARG_EXPECTED_LIB_VARS)
        list(APPEND _ALL_OE_LIB_VARS ${ARG_EXPECTED_LIB_VARS})
        list(REMOVE_DUPLICATES _ALL_OE_LIB_VARS)
    endif()

    set(_EXPECTED_LIBS "")
    if(OpenEye_LIBRARY_TYPE STREQUAL "SHARED")
        foreach(_LIB_VAR ${_ALL_OE_LIB_VARS})
            if(${_LIB_VAR})
                get_filename_component(_LIB_NAME "${${_LIB_VAR}}" NAME)
                if(_EXPECTED_LIBS)
                    set(_EXPECTED_LIBS "${_EXPECTED_LIBS}, \"${_LIB_NAME}\"")
                else()
                    set(_EXPECTED_LIBS "\"${_LIB_NAME}\"")
                endif()
            endif()
        endforeach()
    endif()

    set(_BUILD_INFO_CONTENT "# Auto-generated build information - do not edit
OPENEYE_BUILD_VERSION = \"${_BUILD_INFO_VERSION}\"
OPENEYE_LIBRARY_TYPE = \"${OpenEye_LIBRARY_TYPE}\"
OPENEYE_EXPECTED_LIBS = [${_EXPECTED_LIBS}]
")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/_build_info.py "${_BUILD_INFO_CONTENT}")
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/_build_info.py
        DESTINATION ${ARG_NAME}
    )

    # Also copy _build_info.py to the output directory for editable installs
    add_custom_command(TARGET ${_TARGET_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_BINARY_DIR}/_build_info.py
            ${ARG_PYTHON_OUTPUT_DIR}/_build_info.py
        COMMENT "Copying _build_info.py to ${ARG_PYTHON_OUTPUT_DIR}"
    )

endfunction()
