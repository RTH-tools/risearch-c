cmake_minimum_required(VERSION 3.20)
project(risearch1 LANGUAGES C CXX)

# Development and testing build for RIsearch1 only. The shipped wheel is still
# built by src/Makefile via setup.py, which this does not touch.

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# Path RIsearch searches for on-disk energy tables. RIsearch1 has its matrices
# compiled in and never reads this, but the source requires the macro to exist.
set(RISEARCH_MATPATH "${CMAKE_CURRENT_SOURCE_DIR}/../tables" CACHE STRING
    "Directory holding energy matrix files")

set(RISEARCH1_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src)
file(GLOB RISEARCH1_SOURCES CONFIGURE_DEPENDS ${RISEARCH1_SRC}/*.c)

# The RIsearch1 translation units as a library, so tests can call into them.
# main() is renamed rather than excluded: it lives in risearch.c alongside the
# functions under test, and the test binary supplies its own main(). The rename
# is PRIVATE so it applies only to these sources, not to anything linking them.
add_library(risearch1_core STATIC ${RISEARCH1_SOURCES})
target_compile_definitions(risearch1_core
    PUBLIC  MATPATH="${RISEARCH_MATPATH}" RISVERSION=1
    PRIVATE main=risearch1_main
)
target_include_directories(risearch1_core PUBLIC ${RISEARCH1_SRC})
target_link_libraries(risearch1_core PUBLIC m)

# The CLI, compiled from the same sources but without the main() rename, so no
# extra shim file has to sit in RIsearch1/src where the Makefile would glob it.
# Output lands in the build tree; the Makefile's RIsearch1/bin/RIsearch is left
# alone.
add_executable(RIsearch ${RISEARCH1_SOURCES})
target_compile_definitions(RIsearch PRIVATE
    MATPATH="${RISEARCH_MATPATH}" RISVERSION=1)
target_link_libraries(RIsearch PRIVATE m)

add_subdirectory(tests)
