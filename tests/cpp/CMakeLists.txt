include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.15.2
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(googletest)

add_executable(oemmpa_tests
    test_errors.cpp
    test_molecule_record.cpp
    test_load_report.cpp
    test_oemmpa.cpp
    test_data_objects.cpp
    test_pair_scoring.cpp
    test_statistics.cpp
    test_query_options.cpp
    test_fragmentation_strategy.cpp
    test_fragmenter.cpp
    test_fragment_smiles.cpp
    test_ring_perception.cpp
    test_stereo_enumeration.cpp
    test_environment_fingerprint.cpp
    test_mcs_common.cpp
    test_memory_index.cpp
    test_series.cpp
    test_analyzer.cpp
    test_transform_application.cpp
    test_wizepairz.cpp
)

if(OEMMPA_HAS_DUCKDB)
    target_sources(oemmpa_tests PRIVATE test_duckdb_store.cpp)
endif()

target_link_libraries(oemmpa_tests
    PRIVATE
        oemmpa
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(oemmpa_tests)

add_executable(oemmpa_prefilter_dump prefilter_parity_dump.cpp)
target_link_libraries(oemmpa_prefilter_dump PRIVATE oemmpa)

add_executable(oemmpa_prefilter_dump_disabled prefilter_parity_dump.cpp)
target_link_libraries(oemmpa_prefilter_dump_disabled PRIVATE oemmpa_no_prefilter)
# The dump translation unit holds no prefilter code, so it does not need the
# disabling macro -- only the same header configuration as the library it links.
target_include_directories(oemmpa_prefilter_dump_disabled PRIVATE
    $<TARGET_PROPERTY:oemmpa,INTERFACE_INCLUDE_DIRECTORIES>)
target_compile_definitions(oemmpa_prefilter_dump_disabled PRIVATE
    $<TARGET_PROPERTY:oemmpa,INTERFACE_COMPILE_DEFINITIONS>)

add_test(NAME WizePairZPrefilterParity
    COMMAND ${CMAKE_COMMAND}
        -DENABLED_EXE=$<TARGET_FILE:oemmpa_prefilter_dump>
        -DDISABLED_EXE=$<TARGET_FILE:oemmpa_prefilter_dump_disabled>
        -DFIXTURE=${CMAKE_CURRENT_SOURCE_DIR}/../data/surechembl_mmp_fixture.smi
        -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/prefilter_parity.cmake
)
