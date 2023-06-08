from skbuild import setup  # type: ignore
import sys

# This seems to be needed to ensure Python uses the correct gcc. Without this
# You may get a linker warning, for example:
#       ld: warning: dylib (/usr/local/lib/librssringoccs.so) was built for
#       newer macOS version (10.15) than being linked (10.9)
# librssringoccs is built using a shell script and makes calls to gcc, whereas
# this file is using Python's distutils to handle the compiling. It seems
# Python may use the wrong compiler, causing this error. Setting the following
# CFLAG fixed the issue.

# #   We only need this fix for macOS, so check what operating system is used.
# if (platform.system() == "Darwin"):
#     os.environ["CFLAGS"] = "-mmacosx-version-min=%s" % platform.mac_ver()[0]
#
# srclist = []
#
# for file in os.listdir("rss_ringoccs/crssringoccs"):
#
#     # Only add .c files.
#     if (file[-1] == "c"):
#         srclist.append("rss_ringoccs/crssringoccs/%s" % file)

setup(
    name="rss_ringoccs",
    version="1.3",
    description="C Tools for rss_ringoccs",
    author="Ryan Maguire",
    packages=[
        "rss_ringoccs",
        "rss_ringoccs.calibration",
        "rss_ringoccs.occgeo",
        "rss_ringoccs.rsr_reader",
        "rss_ringoccs.scatter",
        "rss_ringoccs.tools",
    ],
    url="https://github.com/NASA-Planetary-Science/rss_ringoccs",
    project_urls={
        "Documentation": "https://rss-ringoccs.readthedocs.io/en/master/",
        "Bug Reports": "https://github.com/NASA-Planetary-Science/rss_ringoccs/issues",
        "Source": "https://github.com/NASA-Planetary-Science/rss_ringoccs",
    },
    install_requires=[
        "cmake",
        "numpy",
        "scipy",
        "spiceypy",
        "matplotlib",
        "mayavi",
        "pandas",
        "PyMieScatt",
    ],
    cmake_install_dir="rss_ringoccs",
    cmake_args=[
        "-DBUILD_EXAMPLES:BOOL=OFF",
        "-DBUILD_TESTS:BOOL=OFF",
        "-DPython_VERSION={}.{}.{}".format(
            sys.version_info.major, sys.version_info.minor, sys.version_info.micro
        ),
    ],
)
