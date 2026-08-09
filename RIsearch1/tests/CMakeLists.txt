include(FetchContent)

FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        v1.15.2
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Every test and benchmark under tests/, at any depth, so a new file needs no
# CMake change. Same pattern the library uses for src/.
file(GLOB_RECURSE RISEARCH1_TEST_SOURCES CONFIGURE_DEPENDS
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_executable(risearch1_tests ${RISEARCH1_TEST_SOURCES})
target_link_libraries(risearch1_tests PRIVATE risearch1_core GTest::gtest_main)
# The benchmark includes linspace.h directly, so this target compiles the same
# headers the library does and needs the same scaffolding flag.
target_compile_options(risearch1_tests PRIVATE -fpermissive ${RISEARCH1_WARNINGS})
# risearch1_core is built with LTO, so its archive holds IR rather than finished
# objects; this target has to be able to consume that.
risearch1_optimise(risearch1_tests)
target_include_directories(risearch1_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/system)
target_compile_definitions(risearch1_tests PRIVATE
        RISEARCH_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/system/data")