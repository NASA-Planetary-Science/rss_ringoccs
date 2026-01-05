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

# Python is swapping distutils with setuptools.
try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

import os
import platform
import subprocess

# Numpy is needed for the get_include function which tells us where various
# header files, like numpy/arrayobject.h, live.
import numpy

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

add_directory("auxiliary")
add_directory("common")
add_directory("diffraction_correction")
add_directory("extract_csv_data")
add_directory("get_merged_csv_data")
add_directory("get_uranus_data")
add_directory("py_csv_obj")
srclist.append("rss_ringoccs/crssringoccs/crssringoccs.c")

subprocess.call(["make", "-j", "BUILD_STATIC=1"])

setup(
    name = "rss_ringoccs",
    version = "1.3",
    description = "C Tools for rss_ringoccs",
    author = "Ryan Maguire",
    ext_modules = [
        Extension(
            "crssringoccs",
            srclist,
            include_dirs=[
                numpy.get_include(),
                "./",
                "../"
            ],
            extra_objects = [
                "./librssringoccs.a",
                "./libtmpl/libtmpl.a"
            ]
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
