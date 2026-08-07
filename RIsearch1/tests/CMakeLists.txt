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
        test_string_utils.cpp
)
target_link_libraries(risearch1_tests PRIVATE risearch1_core GTest::gtest_main)
