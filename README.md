# Table of Contents
1. [rss_ringoccs](#rss_ringoccs)
    1. [Introduction](#introduction)
    2. [Building](#building)
        1. [OpenMP Support](#openmp-Support)
        2. [Obtaining rss_ringoccs](#obtaining-rss_ringoccs)
        3. [Building (GNU / Linux)](#building-gnu--linux)
        4. [Building (FreeBSD)](#building-freebsd)
        5. [Building (macOS)](#building-macos)
            1. [OpenMP using gcc](#openmp-using-gcc)
            2. [OpenMP using Apple's clang](#openmp-using-apples-clang)
            3. [Building](#building-1)

# rss_ringoccs
rss_ringoccs is a suite of open-source C and Python-based analysis tools for
Cassini Radio Science (RSS) ring occultations. It was developed by Team Cassini
at Wellesley College (Sophia Flury, Jolene Fong, Ryan Maguire, and Glenn
Steranka) under the direction of Richard French, Cassini RSS Team Leader, with
funding provided by the NASA/JPL Cassini project.
## Introduction
The Cassini Radio Science Subsystem (RSS) was used during the Cassini orbital
tour of Saturn to observe a superb series of ring occultations that resulted in
high-resolution, high-SNR radial profiles of Saturn's rings at three radio
wavelengths: 13 cm (S band), 3.6 cm (X band), and 0.9 cm (Ka band). Radial
optical depth profiles of the rings at 1- and 10-km resolution produced by the
Cassini RSS team, using state of the art signal processing techniques to remove
diffraction effects, are available on NASA's Planetary Data System (PDS). These
archived products are likely to be quite adequate for many ring scientists, but
for those who wish to generate their own diffraction-reconstructed ring profiles
from Cassini RSS observations, we offer rss_ringoccs: a suite of Python-based
analysis tools for Cassini Radio Science (RSS) ring occultations.

The purpose of rss_ringoccs is to enable scientists to produce "on demand"
radial optical depth profiles of Saturn's rings from the raw RSS data, without
requiring a deep familiarity with the complex processing steps involved in
calibrating the data and correcting for the effects of diffraction. The code
and algorithms are extensively documented, providing a starting point for users
who wish to test, refine, or optimize the straightforward methods we have
employed. Our emphasis has been on clarity, sometimes at the expense of
programming efficiency and execution time. rss_ringoccs does an excellent job
of reproducing existing 1 km-resolution RSS processed ring occultation data
already present on NASA's PDS Ring-Moons Node, but we make no claim to having
achieved the state-of-the-art in every respect. We encourage users to augment
our algorithms and to report on those improvements, so that they can be
incorporated in future editions of rss_ringoccs.

## Building

`rss_ringoccs` has a few dependencies:

1. A `C` compiler (`gcc`, `clang`, `tcc`, `pcc`, `icx`, and `MSVC` work fine).
2. `python3` (we have tested versions `3.7` to `3.14`).
3. `libtmpl`
4. `setuptools` (or `distutils` if you are using `python 3.10` or older).
5. `numpy`
6. `scipy`
7. `spiceypy`
8. `CSPICE`
9. `GNU Make` (Unix-like) or `CMake` (Windows or Unix-like).
10. `git` (to download the source code).

`libtmpl` is provided as a submodule, and the Python dependencies can be
found in the `requirements.txt` file.

`CPSICE` can be obtained through the
[NAIF website](https://naif.jpl.nasa.gov/naif/toolkit_C.html).
1.  Make sure the CSPICE library files are in your path when building.
    For non-Windows users, this can be achieved by placing
    the files `cspice.a` and `csupport.a` in `/usr/local/lib`.
    For Windows users this can be done by adding the extracted `cspice\lib`
    directory to your `PATH`.
2.  `librssringoccs` provides a small header file with declarations for
    all of the CSPICE functions that are needed
    (`rss_ringoccs` only needs a few).
    You do not need to copy the CSPICE header files to `/usr/local/include`
    or any similar directory.

### OpenMP Support

It is **HIGHLY** recommended that you compile `rss_ringoccs` with OpenMP
support. The inner for-loops in the processing steps can be parallelized,
resulting in a significant speed boost (about 30x on a 32-core CPU).

Parallelizing does require more memory due to thread safety issues
(`n` cores requires `n` times the memory allocated for arrays so that each
thread has its own private data). This increase is relatively small, the
highest resolution reconstructions attempted so far required about 2GB of
memory with 32 threads all firing at once.

Because of this, if you have limited memory (say, less than 8GB of RAM)
and intend to perform high resolution processing, then you should **not**
enable OpenMP.

Note that OpenMP is not a required dependency and both `librssringoccs` and
`libtmpl` will compile with or without OpenMP enabled.
Instructions for installing with OpenMP support are provided for each
supported platform below.

### Obtaining rss_ringoccs

To obtain `rss_ringoccs`, clone the repository:

```bash
git clone --recursive http://github.com/NASA-Planetary-Science/rss_ringoccs.git
cd rss_ringoccs/
```

This will also clone `libtmpl`, which is the primary `C` depenedency.

### Building (GNU / Linux)

On Debian / Ubuntu-based operating systems, you can install the needed tools
using:

```bash
sudo apt install gcc python3 make git
```

If you prefer to use `CMake`, simply do:

```bash
sudo apt install gcc python3 cmake git
```

Most other GNU / Linux systems have these packages readily available.
Use your package manager to obtain them (`pacman`, `yum`, etc.).

If you are using `gcc`, then OpenMP support is already included.
If you are using LLVM's `clang` on Debian GNU / Linux (or similar), install via:

```bash
sudo apt install libomp-dev
```

Similar installation instructions exist for other distributions.
If you are using a compiler other than `gcc` or `clang`, consult the compilers
manual to see if it supports OpenMP and the `-fopenmp` flag.

To build, run the following:

```bash
export USE_OPENMP=1

python3 -m venv .venv
source .venv/bin/activate

python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
python3 -m pip install .
```

If you do not want OpenMP support enabled, omit the `export USE_OPENMP=1` line.

### Building (FreeBSD)

The build instructions for FreeBSD are almost identical to GNU / Linux systems.
On FreeBSD (and other BSDs) you'll need to use `gmake`
(the Makefile uses `GNU Make` features, the default FreeBSD make will not work)
or `CMake`. The build system does **not** support BSD `make`.
Install the required packages with:

```bash
sudo pkg install gcc python3 gmake git
```

Or

```bash
sudo pkg install gcc python3 cmake git
```

To enable OpenMP support, either install `gcc` or a recent version of LLVM:

```bash
sudo pkg install gcc
export CC=gcc
```

or

```bash
sudo pkg install llvm21
```

You may then build `rss_ringoccs` as follows:

```bash
export USE_OPENMP=1

python3 -m venv .venv
source .venv/bin/activate

python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
python3 -m pip install .
```

If you do not want OpenMP support enabled, omit the `export USE_OPENMP=1` line.

### Building (macOS)

On `macOS`, install developer tools with:

```bash
xcode-select --install
```

Enabling OpenMP is definitely more challenging on macOS then the other
supported operating systems, but it can be done using Homebrew.

#### OpenMP using GCC

Perhaps the easiest method is to simply use a different compiler.
Apple's `C` compiler does not support OpenMP out of the box.

1. Install [Homebrew](https://brew.sh/).
2. Install `CMake` and `gcc` using homebrew:
```bash
    brew install cmake gcc
```
3. Set your `CC` environment variable to the `homebrew` version of `gcc`:
```bash
    export CC=$(brew --prefix)/bin/gcc-15
```

Replace this file path with whichever version you installed.

**NOTE:** Do *not* use `/use/bin/gcc`. On macOS this is still Apple's
version of `clang`. OpenMP support will not compile if you use this version.

#### OpenMP using Apple's clang

It is possible to use Apple's `C` compiler with OpenMP, but you'll need to
provide `libomp` yourself. To do this, try the following.

1. Install [Homebrew](https://brew.sh/).
2. Install `CMake` and `libomp` using homebrew:
```bash
    brew install cmake libomp
```
3. Set your `C` environment flags to include `libomp` in their search path.
4. Edit `setup.py` to include the `-lomp` linker flag.

This route is much more of a pain, especially since you need to make sure
either `CMake` or `GMake` are able to search for `libomp`. We do not recommend
this method.

#### Building

You may then build `rss_ringoccs` as follows:

```bash
export USE_OPENMP=1

python3 -m venv .venv
source .venv/bin/activate

python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
python3 -m pip install .
```

If you do not want OpenMP support enabled, omit the `export USE_OPENMP=1` line.

### Building (Windows)

Lastly, there is some experimental Windows support using `CMake` and `MSVC`
(other compilers on Windows should work, but have not been tested).
Visit the following links to obtain the necessary build tools:

1. [CMake](https://cmake.org/download/)
2. [MSVC](https://visualstudio.microsoft.com/vs/features/cplusplus/)
3. [Python](https://www.python.org/downloads/windows/)
4. [Git for Windows](https://git-scm.com/install/windows)

**Note:**
There are many ways to install Python on Windows, including through
Microsoft's app store. None of these alternative installations have been
tested, nor will they be supported by `rss_ringoccs`. Please use the standard
version from the official Python website.

#### OpenMP (Windows)

Microsoft's `MSVC` has support for OpenMP.
If you are using `gcc` via MinGW or something similar, OpenMP support is
also available automatically. With LLVM's `clang` you need to make sure
`libomp` is installed and in your `PATH`.
See
[https://clang.llvm.org/docs/OpenMPSupport.html](https://clang.llvm.org/docs/OpenMPSupport.html)
for details.

#### Compiling (Windows)

For Windows the build instructions are slightly different.

```bash
set USE_OPENMP=1
py -m venv .venv
.venv\Scripts\activate.bat

py -m pip install --upgrade pip
py -m pip install -r requirements.txt
py -m pip install .
```

If you do not have OpenMP support, do not add the `set USE_OPENMP=1`
line at the top.

## Updating

This project updates regularly, as does `libtmpl`. To rebuild, follow the
cleaning instruction below, and then pull the latest changes with:

```bash
git pull --recurse-submodules
```

You may then re-build the latest changes using the previous commands.

## Cleaning / Uninstalling

### Not Windows

Uninstalling `rss_ringoccs` can be done using pip:

```bash
python3 -m pip uninstall rss_ringoccs
make clean
```

Ff you created a virtual environment using `python3 -m venv .venv`,
remove this via:

```bash
rm -rf .venv
```

### Windows

Type the following in the command prompt.

```
rmdir /s build
rmdir /s rss_ringoccs.egg-info
rmdir /s .venv
```

## Release notes
`rss_ringoccs` has changed dramatically over the years from a project that
was primarily written in `Python`, to one that is now mostly written in `C`.
Release notes are contained in
https://github.com/NASA-Planetary-Science/rss_ringoccs/blob/master/ReleaseNotes.md

## Batch Data Processing
Once `rss_ringoccs` has been installed and the necessary data and
ephemeris/geometry files have been downloaded to local storage, as describd in
the `rss_ringoccs User Guide`, users should first test `rss_ringoccs` using the
documented example scripts. Once all is well, users will be able to process a
set of Cassini RSS ring observations in batch mode. To simplify and expedite the
use of `rss_ringoccs`, we provide a Python script that performs the end-to-end
pipeline for a list of files contained in a reference ASCII text file. The
default list is the 1 kHz Cassini RSR files prior to the USO failure, which can
be found in the  `./tables/` directory. This batch script implementation of the
pipeline is located in the `./pipeline/` directory. We suggest running the batch
script using the `yes` command as shown here:
```cd rss_ringoccs_master/pipeline
yes | python e2e_batch.py
```
The `rss_ringoccs User Guide` includes several additional examples of end-to-end
processing scripts, as well as instructions to enable users to construct their
own batch end-to-end scripts.

## Help
If you have trouble with installation or execution of the rss_ringoccs package,
we encourage you to post a issue to
https://github.com/NASA-Planetary-Science/rss_ringoccs/issues. We will attempt
to respond promptly, and ther users will benefit. Alternatively, you can write
email directly to Richard French: rfrench_at_wellesley.edu.

## Citing rss_ringoccs
If you use rss_ringoccs as the basis of a publication, please consider
citing rss_ringoccs using the DOI:10.5281/zenodo.2548947

## Version History

| Version     | Date                  |
|-------------|-----------------------|
| 1.3-beta    | January 12, 2021      |
| 1.2.1       | Septmber 7, 2019      |
| 1.2         | July 3, 2019          |
| 1.1         | February 5, 2019      |
| 1.0         | September 30, 2018    |
| Pre-Release | April 22, 2018        |

## Acknowledgements
This work was supported by the NASA/JPL Cassini mission. We are especially
grateful  to Linda Spilker and Kathryn Weld for their encouragement and support,
and to  RSS Team Member Essam Marouf for developing the diffraction
reconstruction technique that underlies this work.
