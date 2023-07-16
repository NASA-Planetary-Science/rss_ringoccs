#!/bin/bash
if [ ! "$BASH_VERSION" ] ; then
        echo "Please use BASH to run this script ($0)" 1>&2
        exit 1
fi
web="https://atmos.nmsu.edu/pdsd/archive/data/"
# post-USO failure 1 kHz RSR files with uplink station documentation
input=../tables/rsr_1kHz_files_after_USO_failure_withuplink.txt
while IFS= read -r dirfil
do
localdirfil="../data/$dirfil"
webdirfil=$web$dirfil
if [ ! -e $localdirfil ]
then
echo Current File $localdirfil does not locally exist. Fetching.
if curl --output /dev/null --silent --head --fail "$webdirfil"; then
echo "URL exists: $webdirfil"
curl --create-dirs -o $localdirfil $webdirfil
else
echo "URL does not exist: $webdirfil"
exit 1
fi
else
echo Current File "$localdirfil" already exists.
fi
done < "$input"
