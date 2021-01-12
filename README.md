# rss_ringoccs
The master branch is currently v1.3-beta. rss_ringoccs is unable to run on the
latest releases of Python (3.8 or 3.9) and it is recommended that 3.6 is used.
There is active work in rewriting the project C89/C90 compliant code to allow
for future use and ease of backwards and forward compatibility.

rss_ringoccs is a suite of open-source C and Python-based analysis tools for
Cassini Radio Science (RSS) ring occultations. It was developed by Team Cassini
at Wellesley College (Sophia Flury, Jolene Fong, Ryan Maguire, and Glenn
Steranka) under the direction of Richard French, Cassini RSS Team Leader, with
funding provided by the NASA/JPL Cassini project.

Version 1.2 was offically released on July 1, 2019.

Version 1.2.1 was a development version that incorporated changes between
Versions 1.2 and 1.3.

Version 1.3-beta was officially released on January 12, 2021

## Introduction
The Cassini Radio Science Subsystem (RSS) was used during the Cassini orbital
tour of Saturn to observe a superb series of ring occultations that resulted in
high-resolution, high-SNR radial profiles of Saturn's rings at three radio
wavelengths: 13 cm (S band), 3.6 cm (X band), and 0.9 cm (Ka band). Radial
optical depth profiles of the rings at 1- and 10-km resolution produced by the
Cassini RSS team, using state of the art signal processing techniques to remove
diffraction effects, are available on NASA's Planetary Data System (PDS). These
archived products are likely to be quite adequate for many ring scientists, but
for those who wish to generate their own diffraction-reconstructed ring profiles
from Cassini RSS observations, we offer rss_ringoccs: a suite of Python-based
analysis tools for Cassini Radio Science (RSS) ring occultations.

The purpose of rss_ringoccs is to enable scientists to produce "on demand"
radial optical depth profiles of Saturn's rings from the raw RSS data, without
requiring a deep familiarity with the complex processing steps involved in
calibrating the data and correcting for the effects of diffraction. The code
and algorithms are extensively documented, providing a starting point for users
who wish to test, refine, or optimize the straightforward methods we have
employed. Our emphasis has been on clarity, sometimes at the expense of
programming efficiency and execution time. rss_ringoccs does an excellent job
of reproducing existing 1 km-resolution RSS processed ring occultation data
already present on NASA's PDS Ring-Moons Node, but we make no claim to having
achieved the state-of-the-art in every respect. We encourage users to augment
our algorithms and to report on those improvements, so that they can be
incorporated in future editions of rss_ringoccs.

## Installation and Documentation
Detailed installation instructions and full documentation are contained the
`rss_ringoccs User Guide` at https://github.com/NASA-Planetary-Science/rss_ringoccs/tree/master/docs/rss_ringoccs_User_Guide_V1.3.pdf.

For experienced users, we provide a `rss_ringoccs Quick Start` guide located at
https://github.com/NASA-Planetary-Science/rss_ringoccs/tree/master/docs/rss_ringoccs_Quick_Start_V1.3.pdf

Release notes are contained in
https://github.com/NASA-Planetary-Science/rss_ringoccs/blob/master/ReleaseNotes.md

Source code documentation is found at
https://rss-ringoccs.readthedocs.io/en/master/

## Batch Data Processing
Once `rss_ringoccs` has been installed and the necessary data and
ephemeris/geometry files have been downloaded to local storage, as describd in
the `rss_ringoccs User Guide`, users should first test `rss_ringoccs` using the
documented example scripts. Once all is well, users will be able to process a
set of Cassini RSS ring observations in batch mode. To simplify and expedite the
use of `rss_ringoccs`, we provide a Python script that performs the end-to-end
pipeline for a list of files contained in a reference ASCII text file. The
default list is the 1 kHz Cassini RSR files prior to the USO failure, which can
be found in the  `./tables/` directory. This batch script implementation of the
pipeline is located in the `./pipeline/` directory. We suggest running the batch
script using the `yes` command as shown here:
```cd rss_ringoccs_master/pipeline
yes | python e2e_batch.py 
```
The `rss_ringoccs User Guide` includes several additional examples of end-to-end
processing scripts, as well as instructions to enable users to construct their
own batch end-to-end scripts.

## How to (Get) Help
If you have trouble with installation or execution of the rss_ringoccs package,
we encourage you to post a issue to
https://github.com/NASA-Planetary-Science/rss_ringoccs/issues. We will attempt
to respond promptly, and ther users will benefit. Alternatively, you can write
email directly to Richard French: rfrench_at_wellesley.edu.
## Citing rss_ringoccs
If you use rss_ringoccs as the basis of a publication, please consider 
citing rss_ringoccs using the DOI:10.5281/zenodo.2548947

## Acknowledgements
This work was supported by the NASA/JPL Cassini mission. We are especially
grateful  to Linda Spilker and Kathryn Weld for their encouragement and support,
and to  RSS Team Member Essam Marouf for developing the diffraction
reconstruction technique that underlies this work.
