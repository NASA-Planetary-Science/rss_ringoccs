#!/bin/bash

################################################################################
#                                  LICENSE                                     #
################################################################################
#   This file is part of libtmpl.                                              #
#                                                                              #
#   libtmpl is free software: you can redistribute it and/or modify            #
#   it under the terms of the GNU General Public License as published by       #
#   the Free Software Foundation, either version 3 of the License, or          #
#   (at your option) any later version.                                        #
#                                                                              #
#   libtmpl is distributed in the hope that it will be useful,                 #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#   GNU General Public License for more details.                               #
#                                                                              #
#   You should have received a copy of the GNU General Public License          #
#   along with libtmpl.  If not, see <https://www.gnu.org/licenses/>.          #
################################################################################
#   Author:     Ryan Maguire                                                   #
#   Date:       February 13, 2022                                              #
################################################################################

# This file is for developers and internal use. It has many pedantic compiler
# arguments to ensure strict ANSI compliant C code.

# Choose whatever C compiler you want. Tested with gcc, clang, tcc, and pcc
# on GNU/Linux (Debian, Ubuntu, Fedora, and more) and FreeBSD 12.1.
CC=cc
STDVER="-std=c89"
USEOMP=0
ExtraArgs=""

for arg in "$@"; do
    if [ "$arg" == "" ]; then
        break
    elif [[ "$arg" == *"-cc"* ]]; then
        CC=${arg#*=}
    elif [[ "$arg" == *"-std"* ]]; then
        STDVER=$arg
    elif [ "$arg" == "-omp" ]; then
        USEOMP=1
    else
        ExtraArgs="$ExtraArgs ${arg#*=}"
    fi
done

if [ $USEOMP == 1 ]; then
    ExtraArgs="$ExtraArgs -fopenmp"
fi

# Compiler arguments. All of these are supported by gcc, tcc, pcc, and clang.
CArgs1="-pedantic -Wall -Wextra -Wpedantic -Wmissing-field-initializers"
CArgs2="-Wconversion -Wmissing-prototypes -Wmissing-declarations"
CArgs3="-Winit-self -Wnull-dereference -Wwrite-strings -Wdouble-promotion"
CArgs4="-Wfloat-conversion -Wstrict-prototypes -Wold-style-definition"
CArgs5="-I../ -I/usr/local/include -DNDEBUG -g -fPIC -O3 -flto -c"
CompilerArgs="$STDVER $ExtraArgs $CArgs1 $CArgs2 $CArgs3 $CArgs4 $CArgs5"

# Linking arguments.
# -O3 is optimization level 3.
# -I../ means include the parent directory so libtmpl/ is in the path.
# -flto is link time optimization.
# -lm means link against the standard math library.
# -o means create an output.
# -shared means the output is a shared object, like a library file.
if [ $USEOMP == 1 ]; then
    LinkerArgs="-O3 -flto -fopenmp -shared -o librssringoccs.so -lm -ltmpl"
else
    LinkerArgs="-O3 -flto -shared -o librssringoccs.so -lm -ltmpl"
fi

# There may be left-over .so and .o files from a previous build. Remove those
# to avoid a faulty build.
echo "Clearing older files..."
rm -f *.so *.o *.obj *.lib

echo "Compiling librssringoccs..."
echo "    Compiler:"
echo "        $CC"
echo "    Version:"
echo "        $STDVER"
echo "    Compiler Options:"
echo "        $CArgs1"
echo "        $CArgs2"
echo "        $CArgs3"
echo "        $CArgs4"
echo "        $CArgs5"
echo "    Extra Compiler Arguments:"
echo "        $ExtraArgs"

# Loop over all directories in src/ and compile all .c files.
for dir in ./src/*; do
    echo ""
    echo "    Compiling $dir"
    for filename in $dir/*.c; do
        echo "        Compiling: $filename"
        if !($CC $CompilerArgs $filename); then
            exit 1
        fi
    done
done

echo ""
echo "Building libtmpl Shared Object (.so file)"
if !($CC ./*.o $LinkerArgs); then
    exit 1
fi

echo "Cleaning up..."
rm -f *.o

echo "Done"

