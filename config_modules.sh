echo -e "\nBuilding modules...\n"
python setup.py config build_ext --inplace

rm -rf build/

#   Move the compiled shared objects (.so) files to the rss_ringoccs/ folder.
#   This way you can import them directly into python if the rss_ringoccs
#   package is in your path.
mv *.so ./rss_ringoccs/
