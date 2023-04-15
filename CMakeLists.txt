cmake_minimum_required(VERSION 3.20)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(CompilerOptions)
include(ConanHelpers)
conan_instalL(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})

include(FetchLibTmpl)

if(SKBUILD)
    message(STATUS "The project is built using scikit-build")
endif()

if(DEFINED ENV{GITHUB_ENV})
    message(STATUS "The project is being built via Github CI")

    # Stop Python lookup as soon as a version satisfying version constraints is founded
    set(Python_FIND_STRATEGY "LOCATION")

    # Windows only: Perform lookup via environment variables before inspecting the registry
    set(Python_FIND_REGISTRY "LAST")

    # macOS only: look for a matching Python version in the form of a framework
    set(Python_FIND_FRAMEWORK "LAST")
endif()

option(BUILD_EXAMPLES "Build all examples" TRUE)
option(BUILD_TESTS "Build all tests" TRUE)
option(ENABLE_ASAN "Enable build with Clang's Address Sanitizer" FALSE)

set(RSS_RINGOCCS_PARENT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.." CACHE PATH "" FORCE)
set(RSS_RINGOCCS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "" FORCE)

project(rss_ringoccs)

cmake_policy(SET CMP0069 NEW)
include(CheckIPOSupported)
# Optional IPO. Do not use IPO if it's not supported by compiler.
check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT IPO_CHECK_OUTPUT LANGUAGES C)
if(IPO_SUPPORTED)
    message(
        STATUS
        "Interprocedural Optimization/Link-Time Optimization is supported."
    )
else()
    message(
        WARNING
        "Interprocedural Optimization/Link-Time Optimization is not supported: ${IPO_CHECK_OUTPUT}"
    )
endif()

add_subdirectory("src")

if(BUILD_TESTS)
    add_subdirectory("tests")
endif()
