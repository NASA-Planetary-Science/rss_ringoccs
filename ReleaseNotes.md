# Release Notes #

## rss_ringoccs Version 1.0 ##

Release date: 2018 September 28

## Known Issues and Limitations ##
1. rss_ringoccs implements effective radial resolution as defined in Marouf, Tyler, and Rosen 1986 (GIVE REF) eq. xx, using a Kaiser-Bessel alpha=2.5 window function. In order to match the 1-km resolution diffraction-corrected profiles produced by Essam Marouf on the PDS Ring-Moon Systems Node, we must specify an effective resolution of about 0.75 km. The source of this discrepancy appears to be a difference in definitions of effective resolution, but is at this time unresolved.

Workaroud: In order to produce the best match to the RSS diffraction-corrected ring profiles on the PDS, specify in rss_ringoccs a desired resolution 0.75 times that given in the PDS files.

2. Power and frequency calibration GUIs give the following error message under some versions of Python on MacOS systems:

-[NSApplication _setup:]: unrecognized selector sent to instance

*** Terminating app due to uncaught exception 'NSInvalidArgumentException', reason: '-[NSApplication _setup:]: unrecognized selector sent to instance'

Workaround: Use Linux operating system, and post an Issue on the Github page for rss_ringoccs

## Planned Augmentations ##

1. Include calculation of threshold optical depth

2. Detailed tutorials to demonstrate the use of the software

3. Runnable scripts to perform push-button diffraction correction, starting either from raw RSS files or from Essam Marouf's (or our) PDS-style geometry, calibration, and diffraction-limited profiles, at any desired resolution (consistent with the sampling theorem and justified by the SNR), for the full set of RSS occultations at S, X, and Ka-band up to the point of USO failure.

4. Summary plots of elevation vs time of the Cassini spacecraft from each DSN, for each occultation

5. Earth views of each occultation

6. Ability to produce summary PDF files for each occultation, in same format as Essam Marouf's

7. Data catalog query - we will work with the PDS to ensure that our recently-submitted RSS ring occultation observation data catalog is compliant with current PDS search capabilities

8. Implement speed improvements using multiprocessing, Cython, and/or C versions of time-consuming routines.

9. Workarounds that remove or significantly reduce the need for GUI components to optimize and streamline the pipeline.
