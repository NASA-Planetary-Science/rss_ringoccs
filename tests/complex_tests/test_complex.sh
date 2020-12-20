
CC=gcc
CompilerArgs1="-O3 -std=c99 -pedantic -pedantic-errors -Wall -Wextra"
CompilerArgs2="-Wpedantic -Wmisleading-indentation -Wmissing-field-initializers"
CompilerArgs3="-Wmissing-prototypes -Wold-style-definition -Winit-self "
CompilerArgs4="-Wmissing-declarations -Wstrict-prototypes"
CompilerArgs="$CompilerArgs1 $CompilerArgs2 $CompilerArgs3 $CompilerArgs4"

for filename in ./*.c; do
    echo -e "\nTest: $filename"
    $CC $CompilerArgs $filename -o test -lrssringoccs -lrssringoccs_compare
    ./test
    rm -f test
done