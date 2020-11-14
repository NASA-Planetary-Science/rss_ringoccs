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


setup(name='special_functions',
      version='1.3',
      description='Common math and physics functions',
      author='Ryan Maguire',
      ext_modules=[
          Extension('special_functions',
                    ['rss_ringoccs/src/specialfunctionsmodule.c'],
                    include_dirs=[numpy.get_include()],
                    library_dirs=['/usr/local/lib'],
                    libraries=['rssringoccs'])
        ]
     )

#   setup(name='diffrec',
#         version='1.3',
#         description='Diffraction correction and modeling',
#         author='Ryan Maguire',
#         ext_modules=[
#             Extension('diffrec',
#                       ['diffractioncorrectionmodule.c'],
#                       include_dirs=[numpy.get_include(), '../include/',
#                                     '../special_functions/'],
#                       library_dirs=['./'],
#                       libraries=['rssringoccsdiffractioncorrection'])
#           ]
#        )

