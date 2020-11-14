#!/bin/bash

#   Apple aliases gcc to mean clang. If you have installed gcc via homebrew,
#   version 10, then you can use the following to compile with that. Otherwise,
#   this will compile with clang on mac and gcc on linux.

#   osstring=`uname`
#   if [ "$osstring" = "Darwin" ]; then
#   	CC=gcc-10
#   elif [ "$osstring" = "Linux" ]; then
#   	CC=gcc
#   else
#       echo "Operating System not recognized"
#       echo "Only MacOSX and Linux supported"
#       echo "Exiting script"
#       exit 1
#   fi

#   To make life easy, we'll just use gcc to mean whatever the operating
#   system aliased it as. As stated, for MacOS this is clang, and for linux
#   distribution this is actually gcc. However, you can set this to whichever
#   compiler you want. It should be noted that rss_ringoccs has been tested on
#   clang and gcc on both linux and mac systems using strict compiler options to
#   ensure compliance with the ANSI standard. Changing this to a different
#   compiler "should" work, though it hasn't been tested.
CC=gcc

echo -e "\nClearing older files..."
rm -f *.so
rm -f *.o

echo "Compiling rss_ringoccs..."
CompilerArgs1="-std=c89 -ansi -pedantic -pedantic-errors -Wall -Wextra"
CompilerArgs2="-Wpedantic -Wmisleading-indentation"
CompilerArgs3="-Wmissing-prototypes -Wold-style-definition"
CompilerArgs4="-Wstrict-prototypes -I./ -DNDEBUG -g -fPIC -O2 -c"
CompilerArgs="$CompilerArgs1 $CompilerArgs2 $CompilerArgs3 $CompilerArgs4"

echo -e "\n\tCompiler Options:"
echo -e "\t\t$CompilerArgs1"
echo -e "\t\t$CompilerArgs2"
echo -e "\t\t$CompilerArgs3"
echo -e "\t\t$CompilerArgs4"

echo -e "\n\tCompiling rss_ringoccs/src/complex/"
for filename in rss_ringoccs/src/complex/*.c; do
    echo -e "\t\tCompiling: $filename"
    $CC $CompilerArgs $filename
done

echo -e "\n\tCompiling rss_ringoccs/src/math/"
for filename in rss_ringoccs/src/math/*.c; do
    echo -e "\t\tCompiling: $filename"
    $CC $CompilerArgs $filename
done

echo -e "\n\tCompiling rss_ringoccs/src/special_functions/"
for filename in rss_ringoccs/src/special_functions/*.c; do
    echo -e "\t\tCompiling: $filename"
    $CC $CompilerArgs $filename
done

echo -e "\nBuilding rss_ringoccs Shared Object (.so file)"

sharedobjectlist=""
for filename in ./*.o; do
    sharedobjectlist="$sharedobjectlist $filename"
done

$CC $sharedobjectlist -shared -o librssringoccs.so -lm

echo "Moving to /usr/local/lib/librssringoccs.so"
sudo mv librssringoccs.so /usr/local/lib/librssringoccs.so

echo "Copying include/ directory to /usr/local/include/rss_ringoccs/"
includedir="/usr/local/include/rss_ringoccs/"

#   Check if /usr/local/include/rss_ringoccs/ is already a directory. If not
#   then create this via mkdir.
[ ! -d "$includedir" ] && sudo mkdir -p "$includedir/include/"
sudo cp -r rss_ringoccs/include/ "$includedir/include/"

echo "Cleaning up..."
rm -f *.o

echo "Done"

echo -e "\nBuilding modules...\n"
python setup.py config build_ext --inplace

rm -rf build/

#   Move the compiled shared objects (.so) files to the rss_ringoccs/ folder.
#   This way you can import them directly into python if the rss_ringoccs
#   package is in your path.
mv *.so ./rss_ringoccs/
