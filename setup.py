# Python is swapping distutils with setuptools.
try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

import os, platform

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
if (platform.system() == "Darwin"):
    os.environ["CFLAGS"] = "-mmacosx-version-min=%s" % platform.mac_ver()[0]

srclist = []

for file in os.listdir("rss_ringoccs/crssringoccs"):

    # Only add .c files.
    if (file[-1] == "c"):
        srclist.append("rss_ringoccs/crssringoccs/%s" % file)

setup(
    name = "rss_ringoccs",
    version = "1.3",
    description = "C Tools for rss_ringoccs",
    author = "Ryan Maguire",
    ext_modules=[
        Extension(
            "crssringoccs",
            srclist,
            include_dirs=[
                numpy.get_include(),
                "/usr/local/include"
            ],
            library_dirs=["/usr/local/lib"],
            libraries=["tmpl", "rssringoccs"]
        )
    ],
    packages=[
        "rss_ringoccs",
        "rss_ringoccs.calibration",
        "rss_ringoccs.occgeo",
        "rss_ringoccs.rsr_reader",
        "rss_ringoccs.scatter",
        "rss_ringoccs.tools"
    ]
)

