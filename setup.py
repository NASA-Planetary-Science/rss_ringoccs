"""
################################################################################
#                                   LICENSE                                    #
################################################################################
#   This file is part of rss_ringoccs.                                         #
#                                                                              #
#   rss_ringoccs is free software: you can redistribute it and/or              #
#   modify it under the terms of the GNU General Public License as published   #
#   by the Free Software Foundation, either version 3 of the License, or       #
#   (at your option) any later version.                                        #
#                                                                              #
#   rss_ringoccs is distributed in the hope that it will be useful             #
#   but WITHOUT ANY WARRANTY# without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#   GNU General Public License for more details.                               #
#                                                                              #
#   You should have received a copy of the GNU General Public License          #
#   along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.     #
################################################################################
#   Purpose:                                                                   #
#       Compiles all of the C code and build the python package rss_ringoccs.  #
################################################################################
#                                 DEPENDENCIES                                 #
################################################################################
#   1.) numpy                                                                  #
#           Used for working with homogeneous arrays.                          #
#   2.) setuptools                                                             #
#           Replaces distutils since this has been deprecated.                 #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2024/09/09                                                         #
################################################################################
"""

import os
import platform
import subprocess
import shutil

# Python is swapping distutils with setuptools.
try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

# Numpy is needed for the get_include function which tells us where various
# header files, like numpy/arrayobject.h, live.
import numpy

# Check if OpenMP support has been requested.
USE_OPENMP = os.environ.get('USE_OPENMP') == '1'

# This seems to be needed to ensure Python uses the correct gcc. Without this
# You may get a linker warning, for example:
#       ld: warning: dylib (/usr/local/lib/librssringoccs.so) was built for
#       newer macOS version (10.15) than being linked (10.9)
# librssringoccs is built using a shell script and makes calls to gcc, whereas
# this file is using Python's distutils to handle the compiling. It seems
# Python may use the wrong compiler, causing this error. Setting the following
# CFLAG fixed the issue.

#   We only need this fix for macOS, so check what operating system is used.
if platform.system() == "Darwin":
    os.environ["CFLAGS"] = f"-mmacosx-version-min={platform.mac_ver()[0]}"

srclist = []

def add_directory(directory):
    """
        Function for adding all of the C files in a
        given directory to the list of files to be compiled.
    """
    full_path = f'rss_ringoccs/crssringoccs/{directory}/'
    for file in os.listdir(full_path):
        if file[-1] == "c":
            srclist.append(f'{full_path}{file}')

def openmp_args():
    """
        Returns the list of arguments needed by
        the C compiler to enable OpenMP support.
    """
    if not USE_OPENMP:
        return [], []

    # Microsoft's MSVC compiler supports OpenMP 2.0.
    if os.name == "nt":
        return ["/openmp"], []

    # OpenMP support is also possible on macOS (but tricky to set up).
    if platform.system() == "Darwin":

        # The required flags for Apple's clang and GCC differ.
        cc = os.environ.get("CC", "")

        # GCC on macOS does not need -Xpreprocessor. Adding it breaks the
        # linker when building with OpenMP support.
        if "gcc" in cc:
            return ["-fopenmp"], ["-fopenmp"]

        # Apple's Clang on macOS requires -Xpreprocessor. Omiting results
        # in an "unsupported option: -fopenmp" error.
        return ["-Xpreprocessor", "-fopenmp"], ["-Xpreprocessor", "-fopenmp"]

    # On GNU / Linux, GCC has OpenMP support by default, and LLVM's clang
    # has support as well (but you may need to install libomp to use it).
    return ["-fopenmp"], ["-fopenmp"]

def static_library_paths():
    """
        Generates both librssringoccs.a and libtmpl.a and returns
        (rssringoccs.lib and tmpl.lib on Windows), and returns the
        paths to these generated files.
    """

    # Windows users are required to use CMake. GNU, Linux, FreeBSD, and
    # macOS users may also use CMake.
    if shutil.which("cmake"):

        build_dir = "build"
        cmake_args = ["cmake", "-S", ".", "-B", build_dir]

        # Additional flags are needed if OpenMP support is requested.
        if USE_OPENMP:
            cmake_args.append("-DOMP=ON")

        # Run CMake to generate the build files.
        subprocess.run(cmake_args, check = True)

        # Run the build process and compile librssringoccs and libtmpl.
        subprocess.run(
            ["cmake", "--build", build_dir, "--config", "Release", "-j"],
            check=True,
        )

        # Path for Windows users, note the '\' instead of '/'.
        if os.name == "nt":
            return [
                build_dir + "\\Release\\rssringoccs.lib",
                build_dir + "\\libtmpl\\Release\\tmpl.lib"
            ]

        # Path for everyone else.
        return [
            build_dir + "/librssringoccs.a",
            build_dir + "/libtmpl/libtmpl.a"
        ]

    # GNU Make is 'gmake' on FreeBSD. Check for this first.
    if shutil.which("gmake"):
        command = ["gmake", "-j", "BUILD_STATIC=1"]

    # GNU Make is simply 'make' on GNU / Linux and macOS.
    elif shutil.which("make"):
        command = ["make", "-j", "BUILD_STATIC=1"]

    else:
        raise RuntimeError("Neither CMake nor GNU Make found. Install one.")

    # gmake and make have nearly identical build steps. The only difference
    # is the name of the command. This is saved in the 'command' variable.
    if USE_OPENMP:
        command.append("OMP=1")

    subprocess.run(command, check = True)
    return ["./librssringoccs.a", "./libtmpl/libtmpl.a"]

# Add all of the crssringoccs files (the C-Python API wrapper) to the build.
add_directory("auxiliary")
add_directory("common")
add_directory("diffraction_correction")
add_directory("extract_csv_data")
add_directory("get_merged_csv_data")
add_directory("get_uranus_data")
add_directory("py_csv_obj")
srclist.append("rss_ringoccs/crssringoccs/crssringoccs.c")

extra_objects = static_library_paths()
compile_args, linker_args = openmp_args()

setup(
    name = "rss_ringoccs",
    version = "1.3",
    description = "C Tools for rss_ringoccs",
    author = "Ryan Maguire",
    ext_modules = [
        Extension(
            "crssringoccs",
            srclist,
            include_dirs = [
                numpy.get_include(),
                "./",
                "../"
            ],
            extra_objects = extra_objects,
            extra_compile_args = compile_args,
            extra_link_args = linker_args
        )
    ],
    packages = [
        "rss_ringoccs",
        "rss_ringoccs.calibration",
        "rss_ringoccs.occgeo",
        "rss_ringoccs.rsr_reader",
        "rss_ringoccs.scatter",
        "rss_ringoccs.tools"
    ]
)
