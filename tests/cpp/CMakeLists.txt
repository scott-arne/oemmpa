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
    test_query_options.cpp
)

target_link_libraries(oemmpa_tests
    PRIVATE
        oemmpa
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(oemmpa_tests)
