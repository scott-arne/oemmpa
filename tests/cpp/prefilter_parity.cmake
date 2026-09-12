# Runs the WizePairZ pair dump from both builds of the library and requires the
# sorted natural keys to be byte-identical.
set(enabled_dump "${WORK_DIR}/prefilter_enabled.tsv")
set(disabled_dump "${WORK_DIR}/prefilter_disabled.tsv")

execute_process(COMMAND "${ENABLED_EXE}" "${FIXTURE}" "${enabled_dump}"
                RESULT_VARIABLE enabled_status)
if(NOT enabled_status EQUAL 0)
    message(FATAL_ERROR "prefilter-enabled dump failed: ${enabled_status}")
endif()

execute_process(COMMAND "${DISABLED_EXE}" "${FIXTURE}" "${disabled_dump}"
                RESULT_VARIABLE disabled_status)
if(NOT disabled_status EQUAL 0)
    message(FATAL_ERROR "prefilter-disabled dump failed: ${disabled_status}")
endif()

# Without this guard an empty-versus-empty comparison would pass vacuously.
file(SIZE "${enabled_dump}" enabled_size)
if(enabled_size EQUAL 0)
    message(FATAL_ERROR "the fixture produced no pairs; the comparison is vacuous")
endif()

execute_process(
    COMMAND "${CMAKE_COMMAND}" -E compare_files "${enabled_dump}" "${disabled_dump}"
    RESULT_VARIABLE diff_status)
if(NOT diff_status EQUAL 0)
    message(FATAL_ERROR
        "the WizePairZ prefilter changed the pair set:\n"
        "  with prefilter: ${enabled_dump}\n"
        "  without:        ${disabled_dump}")
endif()
