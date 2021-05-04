#!/bin/bash

CC=gcc

echo -e "\nClearing older files..."
rm -f *.so *.o

echo "Compiling rss_ringoccs..."
CompilerArgs1="-std=c89 -ansi -pedantic -pedantic-errors -Wall -Wextra"
CompilerArgs2="-Wpedantic -Wmisleading-indentation -Wmissing-field-initializers"
CompilerArgs3="-Wmissing-prototypes -Wold-style-definition -Winit-self"
CompilerArgs4="-Wmissing-declarations -Wnull-dereference -Wwrite-strings" #-Wconversion -Wdouble-promotion"
CompilerArgs5="-Wstrict-prototypes -I../../ -I./ -DNDEBUG -g -fPIC -O3 -flto -c"
CompilerArgs="$CompilerArgs1 $CompilerArgs2 $CompilerArgs3"
CompilerArgs="$CompilerArgs $CompilerArgs4 $CompilerArgs5"

echo -e "\n\tCompiler Options:"
echo -e "\t\t$CompilerArgs1"
echo -e "\t\t$CompilerArgs2"
echo -e "\t\t$CompilerArgs3"
echo -e "\t\t$CompilerArgs4"
echo -e "\t\t$CompilerArgs5"

for dir in */; do
    echo -e "\n\tCompiling $dir"
    for filename in $dir*.c; do
        echo -e "\t\tCompiling: $filename"
        $CC $CompilerArgs $filename
    done
done

echo -e "\nBuilding rss_ringoccs Shared Object (.so file)"
$CC ./*.o -O3 -flto -shared -o librssringoccs.so -lm

echo "Moving to /usr/local/lib/librssringoccs.so"
sudo mv librssringoccs.so /usr/local/lib/librssringoccs.so

echo "Copying include/ directory to /usr/local/include/rss_ringoccs/"
includedir="/usr/local/include/rss_ringoccs/"

#   Check if /usr/local/include/rss_ringoccs/ is already a directory. If not
#   then create this via mkdir.
[ ! -d "$includedir" ] && sudo mkdir -p "$includedir/include/"
sudo cp ../include/* "$includedir/include/"

echo "Cleaning up..."
rm -f *.o

echo "Done"

echo -e "\nBuilding modules...\n"
python3 setup.py config build_ext --inplace

rm -rf build/

#   Move the compiled shared objects (.so) files to the rss_ringoccs/ folder.
#   This way you can import them directly into python if the rss_ringoccs
#   package is in your path.
mv *.so ./rss_ringoccs/

echo "Done"
