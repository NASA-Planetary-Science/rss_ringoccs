# rss_ringoccs
rss_ringoccs is a suite of Python-based analysis tools for Cassini Radio Science (RSS) ring occultations. It was developed at Wellesley College under the direction of Richard French (Cassini RSS Team Leader) by Research Assistants Glenn Steranka, Ryan Maguire, and Jolene Fong, with funding provided by the NASA/JPL Cassini project.

THIS IS A PRE-RELEASE VERSION. 

The official release of Version 1.0 (and release notes) will be on September 28, 2018.

## Introduction
The Cassini Radio Science Subsystem (RSS) was used during the Cassini orbital tour of Saturn to observe a superb series of ring occultations that resulted in high-resolution, high-SNR radial profiles of Saturn's rings at three radio wavelengths: 13 cm (S band), 3.6 cm (X band), and 0.9 cm (Ka band). Radial optical depth profiles of the rings at 1- and 10-km resolution produced by the Cassini RSS team, using state of the art signal processing techniques to remove diffraction effects, are available on NASA's Planetary Data System (PDS). These archived products are likely to be quite adequate for many ring scientists, but for those who wish to generate their own diffraction-corrected ring profiles from Cassini RSS observations, we offer rss_ringoccs: a suite of Python-based  analysis tools for Cassini Radio Science (RSS) ring occultations.

The purpose of rss_ringoccs is to enable scientists to produce "on demand" radial optical depth profiles of Saturn's rings from the raw RSS data, without requiring a deep familiarity with the complex processing steps involved in calibrating the data and correcting for the effects of diffraction. The code and algorithms are extensively documented, providing a starting point for users who wish to to test, refine, or optimize the straightforward methods we have employed. Our emphasis has been on clarity, sometimes at the expense of programming efficiency and execution time. rss_ringoccs does an excellent job of reproducing existing RSS processed ring occultation data already present on NASA's PDS Ring-Moons Node, but we make no claim to having achieved the state-of-the-art in every respect. We encourage users to augment our algorithms and to report on those improvements, so that they can be  incorporated in future editions of rss_ringoccs. 

## Installation
Detailed installation instructions are contained in https://github.com/NASA-Planetary-Science/rss_ringoccs/blob/master/docs/rss_ringoccs__User_Guide.pdf 

<!--- 
Streamlined installation instructions are provided in https://github.com/NASA-Planetary-;Science/rss_ringoccs/tree/master/docs/Quick_Start_Installation_Guide_for_rss_ringoccs.pdf 
-->

## Documentation
All documentation for the rss_ringoccs package is located in the https://github.com/NASA-Planetary-Science/rss_ringoccs/tree/master/docs
directory. We recommend that you start with rss_ringoccs__User_Guide.pdf.

<!---
Tutorials in the form of Jupyter notebooks are located in https://github.com/NASA-Planetary-Science/rss_ringoccs/tree/master/tutorials
-->

## How to (Get) Help
If you have trouble with installation or execution of the rss_ringoccs package, we encourage you to post a issue to https://github.com/NASA-Planetary-Science/rss_ringoccs/issues. We will attempt to respond promptly, and ther users will benefit. Alternatively, you can write email directly to Richard French: rfrench_at_wellesley.edu.
## Citing rss_ringoccs
If you use rss_ringoccs as the basis of a publication, please consider 
citing rss_ringoccs using the DOI:10.5281/zenodo.1226748

## Known Working Environments
The rss_ringoccs package was develped and tested on several unix-based operating systems, including MacOS 10.13.6 and Linux distributions Ubantu xx.yy, Debian.zz, others....

rss_ringoccs is written in Python 3. We have tested the code under the following Python distributions:

Enthought Deployment Manager -- https://www.enthought.com -- Python 3.5.2
Anaconda --https://anaconda.org/anaconda/python -- Python 3.7.0

## Acknowledgements
This work was supported by the NASA/JPL Cassini mission. We are especially grateful 
to Linda Spilker and Kathryn Weld for their encouragement and support, and to 
RSS Team Member Essam Marouf for developing the diffraction correction technique
that underlies this work.
