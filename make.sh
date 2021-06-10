#!/bin/bash

CC=cc
STDVER="-std=c89"
USEOMP=0

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
        echo "Invalid argument"
        echo "$arg"
        exit 0
    fi
done

if [ $USEOMP == 1 ]; then
    STDVER="$STDVER -fopenmp"
fi

# Name of the created Share Object file (.so).
SONAME="librssringoccs.so"

# Location to store SONAME at the end of building.
SODIR="/usr/local/lib"

if [ $USEOMP == 1 ]; then
    LinkerArgs="-O3 -I/usr/local/include/ -flto -fopenmp -shared -o $SONAME -lm"
else
    LinkerArgs="-O3 -I/usr/local/include/ -flto -shared -o $SONAME -lm"
fi

LinkerArgs="$LinkerArgs -ltmpl"

# Location where the .h files will be stored.
INCLUDE_TARGET=/usr/local/include/rss_ringoccs

if [ $CC == "gcc" ]; then
    CArgs1="$STDVER -pedantic -pedantic-errors -Wall -Wextra -Wpedantic"
    CArgs2="-Wmisleading-indentation -Wmissing-field-initializers -Wconversion"
    CArgs3="-Wmissing-prototypes -Wold-style-definition -Winit-self"
    CArgs4="-Wmissing-declarations -Wnull-dereference -Wwrite-strings"
    CArgs5="-Wdouble-promotion -Wfloat-conversion -Wstrict-prototypes"
    CArgs6="-I/usr/local/include/ -DNDEBUG -g -fPIC -O3 -flto -c"
    CompilerArgs="$CArgs1 $CArgs2 $CArgs3 $CArgs4 $CArgs5 $CArgs6"

#   Clang has different compiler options, so specify those here if using clang.
elif [ $CC == "clang" ]; then
    CArgs1="$STDVER -pedantic -pedantic-errors -Wall -Wextra -Wpedantic"
    CArgs2="-Wmissing-field-initializers -Wconversion -Weverything"
    CArgs3="-Wmissing-prototypes -Wold-style-definition -Winit-self"
    CArgs4="-Wmissing-declarations -Wnull-dereference -Wwrite-strings"
    CArgs5="-Wdouble-promotion -Wfloat-conversion -Wstrict-prototypes"
    CArgs6="-I/usr/local/include/ -DNDEBUG -g -fPIC -O3 -flto -c"
    CompilerArgs="$CArgs1 $CArgs2 $CArgs3 $CArgs4 $CArgs5 $CArgs6"
elif [ $CC == "tcc" ]; then
    CArgs1="$STDVER -pedantic -Wall -Wextra -Wpedantic"
    CArgs2="-Wmisleading-indentation -Wmissing-field-initializers -Wconversion"
    CArgs3="-Wmissing-prototypes -Wold-style-definition -Winit-self"
    CArgs4="-Wmissing-declarations -Wnull-dereference -Wwrite-strings"
    CArgs5="-Wdouble-promotion -Wfloat-conversion -Wstrict-prototypes"
    CArgs6="-I/usr/local/include/ -DNDEBUG -g -fPIC -O3 -flto -c"
    CompilerArgs="$CArgs1 $CArgs2 $CArgs3 $CArgs4 $CArgs5 $CArgs6"
elif [ $CC == "pcc" ]; then
    CArgs1="$STDVER -pedantic -Wall -Wextra -Wpedantic"
    CArgs2="-Wmisleading-indentation -Wmissing-field-initializers -Wconversion"
    CArgs3="-Wmissing-prototypes -Wold-style-definition -Winit-self"
    CArgs4="-Wmissing-declarations -Wnull-dereference -Wwrite-strings"
    CArgs5="-Wdouble-promotion -Wfloat-conversion -Wstrict-prototypes"
    CArgs6="-I/usr/local/include/ -DNDEBUG -g -fPIC -O3 -flto -c"
    CompilerArgs="$CArgs1 $CArgs2 $CArgs3 $CArgs4 $CArgs5 $CArgs6"
elif [ $CC == "cc" ]; then
    CArgs1="$STDVER -pedantic -Wall -Wextra -Wpedantic"
    CArgs2="-Wmisleading-indentation -Wmissing-field-initializers -Wconversion"
    CArgs3="-Wmissing-prototypes -Wold-style-definition -Winit-self"
    CArgs4="-Wmissing-declarations -Wnull-dereference -Wwrite-strings"
    CArgs5="-Wdouble-promotion -Wfloat-conversion -Wstrict-prototypes"
    CArgs6="-I/usr/local/include/ -DNDEBUG -g -fPIC -O3 -flto -c"
    CompilerArgs="$CArgs1 $CArgs2 $CArgs3 $CArgs4 $CArgs5 $CArgs6"
fi

echo -e "\nClearing older files..."
rm -f *.so *.o

if [ -d "$INCLUDE_TARGET" ]; then
    sudo rm -rf "$INCLUDE_TARGET";
fi

if [ -d "$SODIR/$SONAME" ]; then
    sudo rm -f "$SODIR/$SONAME";
fi

echo "Copying include/ directory to /usr/local/include/rss_ringoccs/"
sudo mkdir -p "$INCLUDE_TARGET/include/"
sudo cp ./include/*.h "$INCLUDE_TARGET/include/"

echo "Compiling librssringoccs..."
echo -e "\n\tCompiler Options:"
echo -e "\t\tCompiler: $CC"
echo -e "\t\t$CArgs1"
echo -e "\t\t$CArgs2"
echo -e "\t\t$CArgs3"
echo -e "\t\t$CArgs4"
echo -e "\t\t$CArgs5"
echo -e "\t\t$CArgs6"

#   Loop over all directories in src/ and compile all .c files.
for dir in ./src/*; do
    echo -e "\n\tCompiling $dir"
    for filename in $dir/*.c; do
        echo -e "\t\tCompiling: $filename"
        if !($CC $CompilerArgs $filename); then
            exit 1
        fi
    done
done

echo -e "\nBuilding librssringoccs Shared Object (.so file)"

$CC ./*.o $LinkerArgs

echo "Moving to /usr/local/lib/librssringoccs.so"
sudo mv $SONAME $SODIR

echo "Cleaning up..."
rm -f *.so *.o

LDPATH=/usr/local/lib
if [[ $LD_LIBRARY_PATH == "" ]]; then
    CREATE_NEW_LD_PATH="LD_LIBRARY_PATH=$LDPATH"
    echo -e "\n# Needed for loading librssringoccs." >> ~/.bashrc
    echo "$CREATE_NEW_LD_PATH" >> ~/.bashrc
    echo "export LD_LIBRARY_PATH" >> ~/.bashrc
    source ~/.bashrc
elif [[ $LD_LIBRARY_PATH != *"$LDPATH"* ]]; then
    CREATE_NEW_LD_PATH="LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LDPATH"
    echo -e "\n# Needed for loading librssringoccs." >> ~/.bashrc
    echo "$CREATE_NEW_LD_PATH" >> ~/.bashrc
    echo "export LD_LIBRARY_PATH" >> ~/.bashrc
    source ~/.bashrc
fi


echo -e "\nBuilding modules...\n"
python3 setup.py config build_ext --inplace

rm -rf build/

#   Move the compiled shared objects (.so) files to the rss_ringoccs/ folder.
#   This way you can import them directly into python if the rss_ringoccs
#   package is in your path.
mv *.so ./rss_ringoccs/


