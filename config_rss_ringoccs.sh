#!/bin/bash
if [ ! "$BASH_VERSION" ] ; then
        echo "Please use BASH to run this script ($0)" 1>&2
        exit 1
fi
osstring=`uname`

if [ "$osstring" = "Darwin" ]; then
	BASH_FILE=bash_profile
	fil="Anaconda3-2020.02-MacOSX-x86_64.sh"
elif [ "$osstring" = "Linux" ]; then
	BASH_FILE=bashrc
	fil="Anaconda3-2020.02-Linux-x86_64.sh"
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
if [ -d ~/opt/anaconda3 ]; then
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

# Source the bashrc/bash_profile file.
source .$BASH_FILE

# Check that anaconda3 is included in your path
echo -e ' \t ' "Checking your PATH variable..."
if [[ ":$PATH:" == *"~/opt/anaconda3/bin"* ]]; then
  echo -e ' \t ' "PATH variable is set"
else
	echo -e ' \t ' "Adding ~/opt/anaconda3/bin to your PATH"
	echo 'export PATH="~/opt/anaconda3/bin:$PATH"' >>~/.$BASH_FILE
fi

# Source .bashrc/.bash_profile to update you path.
echo -e ' \t ' "Sourcing .bashrc/.bash_profile"
echo -e ' \t ' "Make sure you are using Bash when running rss_ringoccs"
source .$BASH_FILE

# Update conda
echo -e ' \t ' "Updating conda..."
yes | conda update conda

# Install Spiceypy
echo -e ' \t ' "Installing spiceypy..."
conda config --add channels conda-forge
yes | conda install -c conda-forge spiceypy

# Install PyMieScatt
yes | conda install -c conda-forge pymiescatt

# Source the .bashrc/.bash_profile again after conda updates.
source .$BASH_FILE

# Change back to working directory
cd "$MY_DIRECTORY"

# Install rss_ringoccs.
echo -e ' \t ' "Building librssringoccs..."
cd librssringoccs/
chmod +x config_librssringoccs.sh
./config_librssringoccs.sh

# Compile the new binary files.
echo -e ' \t ' "Building modules..."
cd ../
chmod +x config_modules.sh
./config_modules.sh

echo "Installation complete"
