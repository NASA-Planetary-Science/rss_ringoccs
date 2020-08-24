rm -f *.so
CFLAGS="-lfftw3 -lm" python setup.py config --compiler=gnu99 build_ext --inplace
rm -rf build/
