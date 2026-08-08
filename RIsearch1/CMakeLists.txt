cmake_minimum_required(VERSION 3.20)
project(risearch1 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(RISEARCH_BUILD_TESTS "Build the gtest suite" ON)

# Path RIsearch searches for on-disk energy tables. RIsearch1 has its matrices
# compiled in and never reads this, but the source requires the macro to exist.
set(RISEARCH_MATPATH "${CMAKE_CURRENT_SOURCE_DIR}/../tables" CACHE STRING
    "Directory holding energy matrix files")

set(RISEARCH1_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src)
# Every source under src/, at any depth, so a new subdirectory needs no CMake change.
file(GLOB_RECURSE RISEARCH1_SOURCES CONFIGURE_DEPENDS
     ${RISEARCH1_SRC}/*.c ${RISEARCH1_SRC}/*.cpp)

# Same for the include path: every directory that holds a header goes on it, so
# "Matrix.h" and "math/Matrix.h" both resolve.
file(GLOB_RECURSE RISEARCH1_HEADERS CONFIGURE_DEPENDS
     ${RISEARCH1_SRC}/*.h ${RISEARCH1_SRC}/*.hpp)
set(RISEARCH1_INCLUDES ${RISEARCH1_SRC})
foreach(header ${RISEARCH1_HEADERS})
    get_filename_component(dir ${header} DIRECTORY)
    list(APPEND RISEARCH1_INCLUDES ${dir})
endforeach()
list(REMOVE_DUPLICATES RISEARCH1_INCLUDES)

set(RISEARCH1_DEFINES MATPATH="${RISEARCH_MATPATH}" RISVERSION=1)
# -fpermissive downgrades the void*-to-T* conversions that C allows implicitly and
# C++ does not. It is scaffolding for the C-to-C++ move: each malloc site wants a
# static_cast, and this comes back out once they all have one.
set(RISEARCH1_WARNINGS -Wall -Wextra -Wnon-virtual-dtor -Wduplicated-branches -Wduplicated-cond -pedantic -fpermissive)

# setup.py copies bin/RIsearch into the wheel, so both binaries land where the
# Makefile used to put them rather than in the build tree.
set(RISEARCH1_BIN ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# --- the CLI, optimised: what the wheel ships -------------------------------
add_executable(RIsearch ${RISEARCH1_SOURCES})
target_include_directories(RIsearch PRIVATE ${RISEARCH1_INCLUDES})
target_compile_definitions(RIsearch PRIVATE ${RISEARCH1_DEFINES})
target_compile_options(RIsearch PRIVATE ${RISEARCH1_WARNINGS} -O3)
set_target_properties(RIsearch PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${RISEARCH1_BIN})

# --- the same program with the debug printing compiled in -------------------
add_executable(RIsearch_dbg ${RISEARCH1_SOURCES})
target_include_directories(RIsearch_dbg PRIVATE ${RISEARCH1_INCLUDES})
target_compile_definitions(RIsearch_dbg PRIVATE ${RISEARCH1_DEFINES} DEBUG VERBOSE=2)
target_compile_options(RIsearch_dbg PRIVATE ${RISEARCH1_WARNINGS} -g -O0)
set_target_properties(RIsearch_dbg PROPERTIES
    OUTPUT_NAME RIsearch.dbg
    RUNTIME_OUTPUT_DIRECTORY ${RISEARCH1_BIN})

# --- the same translation units as a library, for the tests -----------------
# main() is renamed so the test binary can link these under gtest's own main().
# The rename is PRIVATE so it applies only to these sources.
add_library(risearch1_core STATIC ${RISEARCH1_SOURCES})
target_include_directories(risearch1_core PUBLIC ${RISEARCH1_INCLUDES})
target_compile_definitions(risearch1_core
    PUBLIC ${RISEARCH1_DEFINES}
    PRIVATE main=risearch1_main)
target_compile_options(risearch1_core PRIVATE ${RISEARCH1_WARNINGS})
target_link_libraries(risearch1_core PUBLIC m)

if(RISEARCH_BUILD_TESTS)
    add_subdirectory(tests)
endif()
