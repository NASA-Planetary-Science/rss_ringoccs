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
cmake_minimum_required(VERSION 3.15)
project(rss_ringoccs C)

option(BUILD_SHARED_LIBS "Build shared library" OFF)
option(NATIVE "Enable native CPU instructions" OFF)
option(OMP "Enable OpenMP support" OFF)

# The Python extension always links librssringoccs and libtmpl statically.
if (SKBUILD)
    set(BUILD_SHARED_LIBS OFF)
endif()

# libtmpl is provided as a submodule. Add this library to the build process.
add_subdirectory(libtmpl)

# librssringoccs also requires CSPICE from NAIF. Note, librssringoccs does
# not require the CSPICE headers, but it does require the compiled library.
set(CMAKE_FIND_LIBRARY_PREFIXES "")
find_library(CSPICE_LIB NAMES cspice libcspice REQUIRED)
find_library(CSUPPORT_LIB NAMES csupport libcsupport REQUIRED)

################################################################################
#                                librssringoccs                                #
################################################################################
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "src/*.c")
add_library(rssringoccs ${SOURCES})
get_filename_component(RSSRINGOCCS_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

if (RSSRINGOCCS_DIR_NAME STREQUAL "rss_ringoccs")
    get_filename_component(
        RSSRINGOCCS_INCLUDE_ROOT ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY
    )
else()
    set(RSSRINGOCCS_INCLUDE_ROOT ${CMAKE_CURRENT_BINARY_DIR}/include_root)
    file(
        COPY include/
        DESTINATION ${RSSRINGOCCS_INCLUDE_ROOT}/rss_ringoccs/include
    )
endif()

# Add libtmpl and rss_ringoccs headers to the include path.
target_include_directories(
    rssringoccs PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${RSSRINGOCCS_INCLUDE_ROOT}>
)

# librssringoccs requires libtmpl and CSPICE. Link these libraries.
target_link_libraries(
    rssringoccs PRIVATE
    tmpl
    ${CSPICE_LIB}
    ${CSUPPORT_LIB}
)

# Add optimization flags and warnings.
if (MSVC)
    target_compile_options(rssringoccs PRIVATE /W4 /O2)
else()
    target_compile_options(
        rssringoccs PRIVATE -Wall -Wextra -Wpedantic -O3 -fPIC
    )
endif()

# libtmpl provides SIMD versions of many routines if the relevant instructions
# are available on the CPU. Enable native CPU instructions if requested.
if (NATIVE)
    if (NOT MSVC)
        include(CheckCCompilerFlag)
        check_c_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
        check_c_compiler_flag("-mcpu=native" COMPILER_SUPPORTS_MCPU_NATIVE)

        if (COMPILER_SUPPORTS_MARCH_NATIVE)
            target_compile_options(rssringoccs PRIVATE -march=native)
        elseif (COMPILER_SUPPORTS_MCPU_NATIVE)
            target_compile_options(rssringoccs PRIVATE -mcpu=native)
        endif()
    endif()
endif()

# Enable OpenMP if requested. This is highly recommended for rss_ringoccs.
# The computationally heavy portions are parallelized when
# OpenMP support is available.
if (OMP)
    find_package(OpenMP)

    if (OpenMP_C_FOUND)
        message(STATUS "rss_ringoccs compiling with OpenMP support.")
        target_link_libraries(rssringoccs PRIVATE OpenMP::OpenMP_C)
    else()
        message(WARNING "OpenMP support requested but not found.")
    endif()
endif()

# Enable link-time optimization (LTO) if supported.
include(CheckIPOSupported)
check_ipo_supported(RESULT IPO_SUPPORTED)

if (IPO_SUPPORTED)
    set_property(TARGET rssringoccs PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# The user may want to install only the C components, and skip the Python
# parts. This can be accomplished by using CMake directly without pip.
# cmake --install will then install the C library only.
install(
    TARGETS rssringoccs
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(DIRECTORY include/ DESTINATION include/rss_ringoccs/include)

################################################################################
#                       crssringoccs (Python extension)                        #
################################################################################
if (SKBUILD)
    find_package(
        Python REQUIRED COMPONENTS Interpreter Development.Module NumPy
    )

    file(
        GLOB_RECURSE CRSSRINGOCCS_SOURCES CONFIGURE_DEPENDS
        "rss_ringoccs/crssringoccs/*.c"
    )

    # crssringoccs makes use of the C-Python API and the Numpy API. It also
    # calls librssringoccs and libtmpl, two C libraries. Link all of these
    # requirements into crssringoccs.
    python_add_library(crssringoccs MODULE WITH_SOABI ${CRSSRINGOCCS_SOURCES})
    target_link_libraries(crssringoccs PRIVATE rssringoccs tmpl Python::NumPy)

    # librssringoccs and libtmpl are built with LTO, so their static libraries
    # contain compiler IR rather than machine code. Link the module with LTO.
    if (IPO_SUPPORTED)
        set_property(
            TARGET crssringoccs PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE
        )
    endif()

    # The C extension module crssringoccs needs to be installed to
    # site packages for the Python part of rss_ringoccs to function properly.
    install(TARGETS crssringoccs LIBRARY DESTINATION . COMPONENT python)
endif()