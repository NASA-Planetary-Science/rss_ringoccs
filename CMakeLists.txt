################################################################################
#                                   LICENSE                                    #
################################################################################
#   This file is part of rss_ringoccs.                                         #
#                                                                              #
#   rss_ringoccs is free software: you can redistribute it and/or modify       #
#   it under the terms of the GNU General Public License as published by       #
#   the Free Software Foundation, either version 3 of the License, or          #
#   (at your option) any later version.                                        #
#                                                                              #
#   rss_ringoccs is distributed in the hope that it will be useful,            #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#   GNU General Public License for more details.                               #
#                                                                              #
#   You should have received a copy of the GNU General Public License          #
#   along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.     #
################################################################################
#   Author:     Ryan Maguire                                                   #
#   Date:       January 23, 2026                                               #
################################################################################
cmake_minimum_required(VERSION 3.14)

project(rss_ringoccs C)

option(BUILD_SHARED_LIBS "Build shared library" OFF)
option(OMP "Enable OpenMP support" OFF)

add_subdirectory(libtmpl)

set(CMAKE_FIND_LIBRARY_PREFIXES "")
find_library(CSPICE_LIB NAMES cspice libcspice REQUIRED)
find_library(CSUPPORT_LIB NAMES csupport libcsupport REQUIRED)

file(GLOB_RECURSE SOURCES "src/*.c")
add_library(rssringoccs ${SOURCES})

target_include_directories(
    rssringoccs PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
    rssringoccs PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
)

# Check for link-time optimization (LTO) support.
include(CheckIPOSupported)
check_ipo_supported(RESULT IPO_SUPPORTED)

# Check for OMP support.
if (OMP)
    find_package(OpenMP)
    if (OpenMP_C_FOUND)
        message(STATUS "rss_ringoccs compiling with OpenMP support.")
    else()
        message(WARNING "OpenMP support requested but not found.")
    endif()
endif()

if (NOT MSVC)

    # Check for -march=native support
    include(CheckCCompilerFlag)
    check_c_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)

    target_compile_options(
        rssringoccs PRIVATE -Wall -Wextra -Wpedantic -O3 -fPIC
    )

    # Add -march=native if supported
    if (COMPILER_SUPPORTS_MARCH_NATIVE)
        target_compile_options(rssringoccs PRIVATE -march=native)
    endif()

    # Add OMP flags if necessary.
    if (OpenMP_C_FOUND)
        if (AppleClang)
            target_compile_options(rssringoccs PRIVATE -Xpreprocessor -fopenmp)
        else()
            target_compile_options(rssringoccs PRIVATE -fopenmp)
        endif()

        target_link_libraries(rssringoccs PRIVATE OpenMP::OpenMP_C)
    endif()
else()
    target_compile_options(rssringoccs PRIVATE /W4 /O2)

    if (OpenMP_C_FOUND)
        target_compile_options(rssringoccs PRIVATE /openmp)
        target_link_libraries(rssringoccs PRIVATE OpenMP::OpenMP_C)
    endif()
endif()

# Enable link-time optimization if supported.
if(IPO_SUPPORTED)
    set_property(TARGET rssringoccs PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

target_link_libraries(
    rssringoccs PRIVATE
    tmpl
    ${CSPICE_LIB}
    ${CSUPPORT_LIB}
)

install(
    TARGETS rssringoccs
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(DIRECTORY include/ DESTINATION include/rss_ringoccs)
