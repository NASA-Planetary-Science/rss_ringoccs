#!/bin/bash
if [ ! "$BASH_VERSION" ] ; then
        echo "Please use BASH to run this script ($0)" 1>&2
        exit 1
fi
osstring=`uname`

if [ "$osstring" = "Darwin" ]; then
	fil="Anaconda3-2019.07-MacOSX-x86_64.sh"
elif [ "$osstring" = "Linux" ]; then
	fil="Anaconda3-2019.07-Linux-x86_64.sh"
else
	echo "Operating System not recognized"
	echo "Only MacOSX and Linux supported"
	echo "Exiting script"
	exit 1
fi

echo "Running rss_ringoccs_config Script:"

web="https://repo.continuum.io/archive/"
webfil="$web$fil"

# Go to your home directory
MY_DIRECTORY=$(pwd)
cd ~
# Check to see if anaconda3 exists on your computer.
if [ -d anaconda3 ]; then
	echo -e ' \t ' "You already have Anaconda on your computer."
	echo -e ' \t ' "Running updates..."
else
	echo -e ' \t ' "Anaconda3 does not currently exist on your computer."
	echo -e ' \t ' "Running setup scripts now..."

	# Use cURL to download Anaconda.
	curl -Ok "$webfil"

	# Run the Anaconda installer.
	chmod +x $fil
	./$fil

	# Remove the shell script from your home directory.
	rm -f $fil
fi

# Source the bashrc file.
source .bashrc

# Check that anaconda3 is included in your path
echo -e ' \t ' "Checking your PATH variable..."
if [[ ":$PATH:" == *"/anaconda3/bin"* ]]; then
  echo -e ' \t ' "PATH variable is set"
else
	echo -e ' \t ' "Adding ~/anaconda3/bin to your PATH"
	echo 'export PATH="~/anaconda3/bin:$PATH"' >>~/.bashrc
fi

# Source .bashrc to update you path.
echo -e ' \t ' "Sourcing .bashrc"
echo -e ' \t ' "Make sure you are using Bash when running rss_ringoccs"
source .bashrc

# Update conda
echo -e ' \t ' "Updating conda..."
yes | conda update conda

# Install Spiceypy
echo -e ' \t ' "Installing spiceypy..."
yes | conda install -c conda-forge spiceypy

# Install PyMieScatt
conda install -c conda-forge pymiescatt

# Source the .bashrc again after conda updates.
source .bashrc

# Change back to working directory
cd "$MY_DIRECTORY"

# Install rss_ringoccs.
cd ./rss_ringoccs/diffrec/

# Remove any old compiled binary files.
rm -f *.so

# Compile the new binary files.
python setup.py config --compiler=gnu99 build_ext --inplace
rm -rf build/

cd "$MY_DIRECTORY"

echo "Installation complete"
