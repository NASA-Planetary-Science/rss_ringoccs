
CC=gcc
CompilerArgs1="-O3 -std=c99 -pedantic -pedantic-errors -Wall -Wextra"
CompilerArgs2="-Wpedantic -Wmisleading-indentation -Wmissing-field-initializers"
CompilerArgs3="-Wmissing-prototypes -Wold-style-definition -Winit-self "
CompilerArgs4="-Wmissing-declarations -Wstrict-prototypes -Wdouble-promotion"
CompilerArgs5="-Wconversion -Wformat-security -Wuninitialized"
CompilerArgs="$CompilerArgs1 $CompilerArgs2 $CompilerArgs3"
CompilerArgs="$CompilerArgs $CompilerArgs4 $CompilerArgs5"
LinkerFlags="-lrssringoccs -lrssringoccs_compare"

echo -e "\nRunning tests with the following compiler arguments:"
echo -e "\t$CompilerArgs1"
echo -e "\t$CompilerArgs2"
echo -e "\t$CompilerArgs3"
echo -e "\t$CompilerArgs4"
echo -e "\t$CompilerArgs5"

for filename in ./*.c; do
    echo -e "\nTest: $filename"
    if [[ $filename == *"erf"* ]] || [[ $filename == *"faddeeva"* ]]; then
        $CC $CompilerArgs $filename -o test $LinkerFlags -lcerf
    else
        $CC $CompilerArgs $filename -o test $LinkerFlags
    fi
    ./test
    rm -f test
done