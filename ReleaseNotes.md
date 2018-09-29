# Release Notes #

## rss_ringoccs Version 1.0 ##

Release date: 2018 September 30

## Known Issues and Limitations ##
1. rss_ringoccs implements effective radial resolution as defined in Marouf, Tyler, and Rosen 1986 (MTR86, Icarus 68, 120-166) eq. 19, using a Kaiser-Bessel alpha=2.5 window function. In contrast, Marouf et al.'s diffraction-reconstructed profiles on the PDS Ring-Moon Systems Node adopt the shortest resolvable wavelength as the
resolution metric. Its inverse is the
highest spatial frequency preserved in the data. The latter is 1 cycle/km for the 1 km
resolution of Marouf's reconstructed profiles. The value corresponds to ~750 m
processing resolution as defined in MTR86. The bandwidth of the lowpass filter in the final stage of the data processing chain determines such frequency and is selected to achieve the desired resolution.

Workaroud: In order to produce the best match to the RSS diffraction-corrected ring profiles on the PDS, specify in rss_ringoccs a desired resolution 0.75 times that given in the PDS files.

2. Power and frequency calibration GUIs give the following error message under some versions of Python on MacOS systems:

-[NSApplication _setup:]: unrecognized selector sent to instance

*** Terminating app due to uncaught exception 'NSInvalidArgumentException', reason: '-[NSApplication _setup:]: unrecognized selector sent to instance'

Workaround: Use Linux operating system, and post an Issue on the Github page for rss_ringoccs

3. Power and frequency calibration GUIs may sometimes not close when users click the "OK" button or the red "X" button.

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

10. Implement low-pass filter as final step of diffraction reconstruction, tuned to give best match to Marouf et al.'s definition of the effective resolution, rather than the MTR86 eq. 19 definition, so that our results will more more directly comparable to Marouf's.
