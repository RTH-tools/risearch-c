include(FetchContent)

FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.15.2
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(risearch1_tests
    test_nucleotide.cpp
    test_operations.cpp
    test_alignment_symbols.cpp
    test_dsm.cpp
    test_fasta.cpp
    system/test_end_to_end.cpp
    system/test_force_start.cpp
        test_string_utils.cpp
)
target_link_libraries(risearch1_tests PRIVATE risearch1_core GTest::gtest_main)
target_include_directories(risearch1_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/system)
target_compile_definitions(risearch1_tests PRIVATE
    RISEARCH_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/system/data")
