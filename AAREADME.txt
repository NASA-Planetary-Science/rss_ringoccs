NASA-Planetary-Science/rss_ringocs/AAREADME.txt

Revisions:
  2018 Mar 22 - rfrench - original version
  
Overview of NASA-Planetary-Science/rss_ringoccs GitHub repository

Purpose:

The purpose of rss_ringoccs is to enable scientists to produce "on demand" high-resolution radial optical depth profiles of Saturn's rings from the raw RSS data, without requiring a deep familiarity with the complex processing steps involved in calibrating the data and correcting for the effects of diffraction. The processing steps and algorithms are extensively documented, and provide a starting point for sophisticated users to test, refine, and optimize the straightforward methods we have employed. Our emphasis has been on clarity, sometimes at the expense of programming efficiency and execution time. Our code does an excellent job of reproducing existing RSS processed ring occultation data already present on NASA's PDS Ring-Moons Node, but we make no claim to having achieved the state-of-the-art in every respect, We encourage users to augment our algorithms and to report on those improvements, so that they can be incorporated in future editions of rss_ringoccs. 

This file contains basic installation instructions for the open-source rss_ringoccs software package to analyze Cassini Radio Science Subsystem (RSS) ring occultation observations. 

This package was produced with support from the Cassini project. 

Developers:
    Richard G. French (Cassini RSS Team Leader) rfrench@wellesley.edu
    Jolene Fong, Research Assistant, Wellesley College
    Ryan Maguire, Research Assistant, Wellesley College
    Glenn Steranka, Research Assistant, Wellesley College

System requirements:
  This software and installation have been tested on the following operating systems:
    Mac OS 10.13.3 
    Linux Ubantu xxxx
    Linux Debian xxxx
  The code is written in Python 3, but is compatible with Python 2.7. It has been tested using the following versions of Python:
    on Mac OS:
        Anaconda (http://www.python.org/) Python 3.6.4
        Enthought(https://store.enthought.com/downloads/) Python 3.5
        Native python 2.7.10 (/usr/bin/python)
  The software requires the installation of SpiceyPy: The NASA JPL NAIF SPICE toolkit wrapper written in Python, available from   
        https://github.com/AndrewAnnex/SpiceyPy
   
The rss_ringoccs package has a strict directory structure that we recommend that users keep intact:
  
  rss_ringoccs_master/
    AAREADME.txt [this file]
    data/         <--- Raw Cassini RSS data ONLY, downloaded from NASA PDS Atmospheres Node 
    docs/         <--- documentation
    kernels/      <--- NASA/JPL/NAIF/SPICE "kernel" files ONLY (spacecraft trajectory files, planetary ephemerides, etc.)
                       downloaded from https://naif.jpl.nasa.gov/pub/naif/
    output/       <--- output files produced by rss_ringoccs software
    src/          <--- source files (Python and shell scripts)

Quick installation guide:

1) Download the rss_ringoccs repository:
    + Navigate to https://github.com/NASA-Planetary-Science/rss_ringoccs
    + Click on the green pulldown button "Clone or Download"
    + Select "Download ZIP" - this will download a file entitled "rss_ringoccs-master.zip"
    
2) Unpack and install the rss_ringoccs repository:
    + Locate directory containing the downloaded zip file - assume that it is ~/Downloads for this example 
    + Open a terminal window and navigate to a directory below which you wish to install the rss_ringocss/ package
      (This should be on a disk with substantial free storage space, if you intend to process a significant volume of
      RSS data, since raw data files are quite large (in excess of 1 GB per wavelength per occultation.) - assume that it is 
      ~/local/ for this example.
    + Use your favorite utility to unpack the zip file. For example:
        cd ~/local/
        unzip ~/Downloads/rss_ringoccs-master.zip
    + This step will create a subdirectory called rss_ringoccs_master/
    
3) Download and install spiceypy from https://github.com/AndrewAnnex/SpiceyPy, following installation instructions 
    in https://github.com/AndrewAnnex/SpiceyPy/blob/master/README.rst
    Test your installation by firing up Python, using this example:
        In[1]: import spiceypy
        In[2]: spiceypy.pi()
        Out[3]: 3.141592653589793
        
Documentation:

rss_ringoccs_master/docs/ contains extensive documentation of the rss_ringocss package, and a QuickStartGuide.txt

Follow the step-by-step instructions in QuickStartGuide.txt to:
  + download selected raw Cassini RSS data required for the demo programs (takes some time)
  + download selected kernel files required for analysis of Cassini RSS data (takes some time)
  + run the non-interactive sample end-to-end processing of a single RSS ring occultation event
  + confirm that your demo results match the expected results
  
Next steps:

  +Read the more extensive documentation End2EndGuide.txt to learn about the sequence of processing steps go from a raw RSS data file        radial optical depth profile of Saturn's rings.
  +Download additional RSS data files and the required kernel files to have a full distribution of all available data (takes considerable    download time and substantial disk space)
  +While you are waiting for the downloads to complete, read the "RSS_Data_Users_Guide.pdf" for additional background on the Cassini RSS      ring occultation experiments.
  
Getting Help:

Feedback is always welcomed! If you having trouble installing rss_occs, please post a report under the "Issues" tab on the rss_ringoccs GitHub page. We will attempt to respond promptly. Or, you can write directly to Richard French (rfrench@wellesley.edu), Cassini Radio Science Team Leader, who led the effort to produce rss_ringoccs.




  
        
    
