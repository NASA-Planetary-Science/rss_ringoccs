# Release Notes #


## rss_ringoccs V1.3 ##

Release date: 2019 August 30

### Changes from V1.2 ###

1. Fix minor typographical errors in user Guide.

2. Clearer explanation of differences between end-to-end batch processing scripts.

3. C code implemented for diffraction reconstruction and special functions.

4. Support restricted to Anaconda Python for compatibility reasons.

5. Enhance post-USO support. This includes addressing the discontinuous frequency offsets and modeling the phase echo introduced by diffraction of the uplink signal by the ring system.

6. Runnable scripts to perform push-button diffraction correction, starting either from raw RSS files or from Essam Marouf's (or our) PDS-style geometry, calibration, and diffraction-limited profiles, at any desired resolution (consistent with the sampling theorem and justified by the SNR), for the full set of RSS occultations at S, X, and Ka-band up to the point of USO failure.

## rss_ringoccs V1.2 ##

Release date: 2019 July 1

### Changes from V1.1 ###

1. Improved sigma-clipping and frequency offset fitting for speed, accuracy, and usability.

2. Adapted code to process Voyager 2 radio occultation data for Uranus' rings.

**Resolved V1.1-1**

1. Introduced push-button scripts starting from either raw RSS files (the `e2e_batch.py` script) or from Essam Marouf's (or our) PDS-style geometry, calibration, and diffraction-limited profiles (the `quick_look.py` script), at any desired resolution (consistent with the sampling theorem and justified by the SNR), for the full set of RSS occultations at S, X, and Ka-band up to the point of USO failure.

**Resolved V1.1-2**

2. Software is now more extensively documented online, in the User's Guide, and within the push-button scripts.

**Resolved V1.1-4**

3. Improved speed of slowest routines in diffraction reconstruction by re-writing in C, wrapping in Python, and utilizing multi-core processing.

**Resolved V1.1-5**

4. Added support for processing many of the post-USO failure RSR files.

**Resolved V1.1-6**

5. Added support for processing the incoherent, or scattered, signal.

### Known Issues and Limitations of V1.2 ###

#### V1.2-1 (carried over from V1.1-1) ####
For the extreme nearly edge-on viewing geometry of Rev133E at X-band, rss_ringoccs gives slightly different results from PDS, traceable to a difference of about 10% in the cubic term of the varaiable psi. The origin of this discrepancy is unknown, but it is not important for any other occultation data sets we have reduced so far, and is relatively minor even for Rev133E at X band.

#### V1.2-2 ####
Some post-USO files contain discontinuous frequency offsets. rss_rings v1.2 does not support the processing of these files.

### Lien list for V1.3 ###

1. Properly proximal egress occultations in Grand Finale orbits

2. Data catalog query - we will work with the PDS to ensure that our recently-submitted RSS ring occultation observation data catalog is compliant with current PDS search capabilities.

3. Improve processing and documenting of the scattered signal. This should include appropriate Doppler footprint contours and correcting the observed frequency drift to account for wavelength-dependents and motion of the spacecraft.


## rss_ringoccs V1.1 ##

Release date: 2019 Feb 1

### Known Issues and Limitations of V1.1 ###

#### V1.1-1 ####
For the extreme nearly edge-on viewing geometry of Rev133E at X-band, rss_ringoccs gives slightly different results from PDS, traceable to a difference of about 10% in the cubic term of the variable psi. The origin of this discrepancy is unknown, but it is not important for any other occultation data sets we have reduced so far, and is relatively minor even for Rev133E at X band.

### Lien list for V1.2 ###

1. Runnable scripts to perform push-button diffraction correction, starting either from raw RSS files or from Essam Marouf's (or our) PDS-style geometry, calibration, and diffraction-limited profiles, at any desired resolution (consistent with the sampling theorem and justified by the SNR), for the full set of RSS occultations at S, X, and Ka-band up to the point of USO failure.

2. More extensive documentation to demonstrate the use of the software.

3. Data catalog query - we will work with the PDS to ensure that our recently-submitted RSS ring occultation observation data catalog is compliant with current PDS search capabilities.

4. Improve speed of slowest routines by using tested multiprocessor code.

5. Explore possibility of processing post-USO failure RSR files.

6. Explore feasibility, level of effort, and value of archiving scattered signal data -- perhaps as a PDART proposal.


## rss_ringoccs V1.0 ##

Release date: 2018 September 30

### Known Issues and Limitations of V1.0 ###
#### V1.0-1 ####
rss_ringoccs implements effective radial resolution as defined in Marouf, Tyler, and Rosen 1986 (MTR86, Icarus 68, 120-166) eq. 19, using a Kaiser-Bessel alpha=2.5 window function. In contrast, Marouf et al.'s diffraction-reconstructed profiles on the PDS Ring-Moon Systems Node adopt the shortest resolvable wavelength as the
resolution metric. Its inverse is the
highest spatial frequency preserved in the data. The latter is 1 cycle/km for the 1 km
resolution of Marouf's reconstructed profiles. The value corresponds to ~750 m
processing resolution as defined in MTR86. The bandwidth of the lowpass filter in the final stage of the data processing chain determines such frequency and is selected to achieve the desired resolution.

Workaroud: In order to produce the best match to the RSS diffraction-reconstructed ring profiles on the PDS, specify in rss_ringoccs a desired resolution 0.75 times that given in the PDS files.

#### V1.0-2 ####
Power and frequency calibration GUIs give the following error message under some versions of Python on MacOS systems:

-[NSApplication _setup:]: unrecognized selector sent to instance

*** Terminating app due to uncaught exception 'NSInvalidArgumentException', reason: '-[NSApplication _setup:]: unrecognized selector sent to instance'

Workaround: Use Linux operating system, and post an Issue on the Github page for rss_ringoccs

#### V1.0-3 ####
Power and frequency calibration GUIs may sometimes not close when users click the "OK" button or the red "X" button.

Workaround: Use Linux operating system, and post an Issue on the Github page for rss_ringoccs
