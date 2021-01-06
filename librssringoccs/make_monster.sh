

CC=gcc

echo -e "\nClearing older files..."
rm -f *.so
rm -f *.o
rm -f monster.c

echo "Compiling rss_ringoccs..."
CompilerArgs1="-std=c89 -ansi -pedantic -pedantic-errors -Wall -Wextra"
CompilerArgs2="-Wpedantic -Wmisleading-indentation -Wmissing-field-initializers"
CompilerArgs3="-Wmissing-prototypes -Wold-style-definition -Winit-self"
CompilerArgs4="-Wmissing-declarations -Wnull-dereference -Wwrite-strings -D__RSS_RINGOCCS_MAKE_MONSTER__"
CompilerArgs5="-Wstrict-prototypes -I../../ -I./ -DNDEBUG -g -fPIC -O3 -c"
CompilerArgs="$CompilerArgs1 $CompilerArgs2 $CompilerArgs3"
CompilerArgs="$CompilerArgs $CompilerArgs4 $CompilerArgs5"

echo -e "\n\tCompiler Options:"
echo -e "\t\t$CompilerArgs1"
echo -e "\t\t$CompilerArgs2"
echo -e "\t\t$CompilerArgs3"
echo -e "\t\t$CompilerArgs4"
echo -e "\t\t$CompilerArgs5"

touch monster.c

echo "#include <stdio.h>" >> monster.c
echo "#include <stdlib.h>" >> monster.c
echo "#include <string.h>" >> monster.c
echo "#include <math.h>" >> monster.c
echo "#include <ctype.h>" >> monster.c
echo "#include <float.h>" >> monster.c

for filename in ../include/*.h; do
    echo "#include \"$filename\"" >> monster.c
done

for dir in */; do
    for filename in $dir*.c; do
        echo "#include \"$filename\"" >> monster.c
    done
done

echo "Compiling"
$CC $CompilerArgs monster.c

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
sudo cp ../include/* "$includedir/include/"

echo "Cleaning up..."
rm -f *.o

echo "Done"
