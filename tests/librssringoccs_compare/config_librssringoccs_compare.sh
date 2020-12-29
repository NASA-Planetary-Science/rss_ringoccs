#!/bin/bash

CC=gcc

echo -e "\nClearing older files..."
rm -f *.so
rm -f *.o

echo "Compiling rss_ringoccs..."
CompilerArgs1="-std=c99 -pedantic -pedantic-errors -Wall -Wextra"
CompilerArgs2="-Wpedantic -Wmisleading-indentation -Wmissing-field-initializers"
CompilerArgs3="-Wmissing-prototypes -Wold-style-definition -Winit-self"
CompilerArgs4="-Wmissing-declarations -Wstrict-prototypes -I./ -O3 -c"
CompilerArgs="$CompilerArgs1 $CompilerArgs2 $CompilerArgs3 $CompilerArgs4"

echo -e "\n\tCompiler Options:"
echo -e "\t\t$CompilerArgs1"
echo -e "\t\t$CompilerArgs2"
echo -e "\t\t$CompilerArgs3"
echo -e "\t\t$CompilerArgs4"

echo -e "\n\tCompiling Test Functions:"
for filename in ./*.c; do
    echo -e "\t\tCompiling: $filename"
    $CC -Wconversion -Wdouble-promotion $CompilerArgs $filename
done

echo -e "\nBuilding Shared Object (.so file)"

sharedobjectlist=""
for filename in ./*.o; do
    sharedobjectlist="$sharedobjectlist $filename"
done

$CC $sharedobjectlist -shared -o librssringoccs_compare.so -lrssringoccs -lm

echo "Moving to /usr/local/lib/librssringoccs_compare.so"
sudo mv librssringoccs_compare.so /usr/local/lib/librssringoccs_compare.so

echo "Copying header file to /usr/local/include/rss_ringoccs/"
includedir="/usr/local/include/rss_ringoccs/"

#   Check if /usr/local/include/rss_ringoccs/ is already a directory. If not
#   then create this via mkdir.
[ ! -d "$includedir" ] && sudo mkdir -p "$includedir/include/"
sudo cp rss_ringoccs_compare_funcs.h "$includedir/include/"

echo "Cleaning up..."
rm -f *.o

echo "Done"
