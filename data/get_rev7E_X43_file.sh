#!/bin/bash
if [ ! "$BASH_VERSION" ] ; then
	echo "Please use BASH to run this script ($0)" 1>&2
	exit 1
fi
web="https://atmos.nmsu.edu/PDS/data/"
dirfil="cors_0105/sroc1_123/rsr/s10sroe2005123_0740nnnx43rd.2a2"
IFS = read -r dirfil
webdirfil=$web$dirfil
if [ ! -e $dirfil ]
then
        echo Current File $dirfil does not locally exist. Fetching.
        curl --create-dirs -o $dirfil $webdirfil
else
        echo Current File "$dirfil" already exists.
fi