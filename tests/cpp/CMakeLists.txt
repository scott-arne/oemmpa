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
    test_desalter.cpp
    test_pair_scoring.cpp
    test_query_options.cpp
    test_fragmentation_strategy.cpp
    test_fragmenter.cpp
    test_environment_fingerprint.cpp
    test_memory_index.cpp
    test_analyzer.cpp
    test_transform_application.cpp
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
