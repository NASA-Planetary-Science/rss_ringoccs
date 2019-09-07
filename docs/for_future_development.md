# Brief overview and ideas for future development of `rss_ringoccs` v1.3

## Accessibility

- The `rss_ringoccs` software is currently hosted on GitHub repository, available by clone or by download
- Data and kernel reference files live on NASA/PDS node server — will need to curl (or wget if supported) to obtain these (very large). Scripts are provided for this purpose.
- **Lien:** make software installable via pip or conda — should help address many of the issues encountered in attempting to use outside of the source directory because this would enable the pipeline to run anywhere. Possible issues to address: location of input files (RSR and kernels), configuration scripts, and 

## The rss_ringoccs pipeline
The source code for the `rss_ringoccs` software is stored in the `./rss_ringoccs/` directory within the top-level directory cloned or downloaded from GitHub.

### RSRReader class
    Contact: J Fong
- Stored in `./rss_ringoccs/rsr_reader/` 
- Files are binary in a NASA PDS format — must be custom-read
- Code is hard-wired to follow NASA Planetary Data System formatting — will parse a header and then read in subsequent data
- **Lien:**
    -  Recover uplink DSN station from header for post-USO occultations
    - Determine nature of post-USO files ending in “5”
    - Investigate files which lack proper coefficients for computing predicted sky frequencies (stored in header but sometimes are NaN, e.g. Rev 82)

### Geometry
    Contact: J Fong
- Stored in `./rss_ringoccs/occgeo/`
- Uses [`spiceypy`](https://spiceypy.readthedocs.io/en/master/), a python version of the [NAIF/SPICE](https://naif.jpl.nasa.gov/naif/toolkit.html) routines, to compute telemetry from NAIF/SPICE reference files called kernel files
- Supports computing ring intercept radius for a specific reference frame (e.g. the F ring)
- Outputs a .TAB and .LBL file. LBL file is a reference file explaining the contents of the .TAB file.
- Not all geometry results are stored in DLP or TAU files, so this file (and object class) is very useful to have
- The freespace regions computed are stored in the LBL file

### Calibration
    Contact: S Flury
- Stored in `./rss_ringoccs/calibration/` 
- Phase steering
    - Calibration instantiates the FreqOffsetFit (`freq_offset_fit.py`) object class
    - FreqOffsetFit calls `calc_freq_offset` — this computes the carrier signal frequency relative to bandwidth center using FFTs and CFTs sampling small windows of data over the ring system and some surrounding freespace
    - FreqOffsetFit calls `calc_f_sky_recon`, a python version of Nicole Rappaport’s PREDICTS fortran code, to compute the spacecraft frequency as seen at the DSN station using the spacecraft telemetry.
    - FreqOffsetFit performs sigma clipping (`create_mask`, `__sigma_clip`, and `__neighbor_check`) to exclude unusable frequencies (B ring obstruction, noise, etc)
    - If post-USO, FreqOffsetFit looks for discontinuities. If some exist, the following step is done for each segment bounded by the discontinuities. This is new and sometimes a bit fickle.
- FreqOffsetFit uses an iterative F-test (`calc_poly_order`) to find the best polynomial order (&le; 9) for fitting the carrier frequency drift. Then, FreqOffsetFit stores the best-fit polynomial of this order.
    - FreqOffsetFit saves a plot of the frequency offset fit results for user evaluation/reference
    - Calibration uses the results stored in the instantiated FreqOffsetFit object to perform “digital steering”, which corrects phase drift in the signal due to changes in the frequency offset (Eqn 2 in MTR86) using the correct_IQ method
    - **Lien:**
        - investigate possible approaches to 2nd order phase steering
        - trouble-shoot and optimize the discontinuity finder
- Power normalization
    - Calibration instantiates the Normalization object class
    - Normalization uses the freespace regions found by Geometry to determine the viable freespace power. 
    - If insufficient freespace power exists, Normalization fits the entire power profile using a first-order polynomial and kicks the user into an interactive mode (superceded by calling the pipeline with “yes | “ to pipe yes into all input, accepting the default fit)
    - Otherwise, Normalization uses its `create_mask` method to limit diffraction pattern effects on the fit.
    - Then, Normalization fits the total freespace power using a low (&le; 3) order polynomial.
    - Normalization saves a plot of the total profile fit with subplots showing the fit in each individual freespace region.
    -  **Lien:** 
        - add support for a cubic spline fit and a piecewise function fit (this has been started but not fully implemented or tested)
        - Calibration saves the frequency offset and power normalization fits a TAB and LBL file.
- **Lien:** full post-USO support

### DiffractionLimitedProfile (DLP for short)
    Contact: J Fong (primary) or S Flury (tau thresh)
- Stored in `./rss_ringoccs/calibration/`
- This takes the calibrated (phase-steered and normalized) data sampled with respect to time and resamples it with respect to ring intercept radius in resample_IQ using results from Geometry.
- DLP computes uncertainties in normal optical depth and signal power by calling calc_tau_thresh to compute the signal-to-noise ratio and threshold optical depth following MTR86 Eqns 23-26.
- After resampling, DiffractionLimitedProfile determines if the occultation is a chord or diametric (if chord, then this splits the chord at rho intercept velocity of zero km/s) and then uses a class method to create a DiffractionLimitedProfile object instance for each profile direction.
- Each DLP class interpolates the Geometry results to resample them to the ring intercept radius of the newly resampled calibrated data. This includes some checks on the interpolation results to prevent negative power or out-of-bounds phase values.
- Each DLP class writes a TAB and LBL file corresponding to the appropriate profile direction.
- **Lien:**
    - Computing upper and lower uncertainty limits for tau thresh and stationary phase following MTR Eqn 27 for inclusion in the DLP TAB or LBL file
    - Proximal orbit support?

### DiffractionCorrection 
    Contact: R Maguire
- Stored in `./rss_ringoccs/diffrec/`
- This is the crowning achievement of the ring occultation radio science experiment: only with the complex signal can diffraction effects be reduced or eliminated. The Fresnel transform requires the signal phase in order to perform this correction, the phase coming from arctan(Q/I).
- The approach to solving the Fresnel transform follows MTR86. This deviates from MTR86 in computing the phase. While MTR86 uses a second order approximation to the phase term ψ (see their Appendix A), DiffractionCorrection allows for higher order approximations of ψ (following MTR86 Appendix A for a fourth-order approximation) and a full solution using Newton-Raphson root-finding (i.e., finding where dψ/dφ = 0 with φ being the stationary phase computed from the calibrated complex signal) to determine the complete expression.
- Unlike the rest of rss_ringoccs, DiffractionCorrection is written in C and wrapped in python to allow for a more hospitable computation time.


### Scatter
    Contact: S Flury
- Stored in `./rss_ringoccs/scatter/` 
- Examines the calibrated (but not resampled or reconstructed) profile
- Uses a STFT with a Hamming window to document the incoherent signal
- Stacks STFT to improve SNR of faint scattered-signal features
- **Lien:** 
    - “straighten”/correct incoherent signal frequency drift using Doppler footprint, Keplerian orbits, and band wavelength so that scattered signal is always aligned with the coherent signal. J Fong has done some preliminary work on the Doppler footprint — she is the contact person for this. See also Chp 4 of F Thomson’s PhD dissertation, Eqn 4.1 and discussion in S4.1.2.
    - investigate the particle size distribution implied by the incoherent signal

## Using the pipeline

- Owing to the division of labor and modular architecture of the software, each step is executed separately. Each step is an object class which, when instantiated, performs the requisite tasks of that step and stores the results as object attributes. Some results are output to TAB and LBL files, which follow the format of the previously published files on the PDS.
- Each step in the pipeline is contingent on the previous steps, requiring a script to instantiate each object in succession. 
- The first four steps (reading, geometry, calibration, resampling to obtain the DLP) need only be run once. DiffractionCorrection contains numerous options, such as window function, Legendre polynomial approximation, radius range, and reconstruction resolution. This step is likely one users will want to run multiple times.
- The information needed by DiffractionCorrection can be assembled from the TAB files output by the previous steps, and a tool--ExtractCSVData--is provided for this purpose.
- Example end-to-end (in the `./examples/` directory)
    - The `e2e_run.py` script serves as an example of the pipeline from reading a raw RSR file to outputting a reconstructed profile
    - The `quick_look_run.py` script serves as an example of the diffraction correction from reading the .TAB files to outputting a reconstructed profile
- Batches (in the `./pipeline/` directory)
    - DLP batch (`dlp_batch.py`) which produces 30 m DLP profiles for all pre-USO and a large portion of post-USO occultations. This only need ever be run once for each file list (one for the pre-USO RSR files and one for the post-USO RSR files). This will allow reconstruction at resolutions from 100 m to 3 km without re-running the initial steps of the pipeline.
    - DiffRec batch (`diffrec_batch.py`) which produces on-demand reconstructed profiles at custom reconstruction resolutions down to 100 m for pre-USO, 500 m for post-USO, for any desired radius range, window function, approximation of ψ.
    - End-to-end batches of various resolutions (`e2e_batch_1km.py`, `e2e_batch_500m.py`, and `e2e_batch_postUSO_1km.py`) - these produce GEO, CAL, DLP, and TAU files corresponding to the processing of raw data contained in RSR files into reconstructed optical depth profiles

## Tables
Tables of reference information are stored in the `./tables/` directory. These include lists of kernel files to download for processing Cassini data, raw RSR files recommended for download and processing, the data catalog, and orbital elements pertaining to gaps in Saturn's rings.

## Science Tools

### Particle Size Distribution
    Contact: S Flury
- Stored in `./rss_ringoccs/tools/dtau_partsize.py`
- Uses Mie scattering physics to predict differential optical depths over a range of particle sizes assuming a power law size distribution with a range of indices
- Outputs predictions to a CSV for user to do analysis
- Contains script for plotting
- **Lien:** 
    - incorporate code to estimate differential optical depth from a set of reconstructed optical depth profiles
    - incorporate particle size estimates from incoherent signal
    - incorporate haze layer test code

### Ring Feature Fitting
    Contact: S Flury
- Stored in `./rss_ringoccs/tools/ringfit.py` 
- Object class which reads in a diffraction-corrected optical-depth profile
- Fits a user-specified feature (given radial limits) with a user-selected function
- Fit is done to 1-exp(-tau), not optical depth or power
- Fit results are stored as attributes
- Error estimates in both radius and time
- Time reported in both year-doy-spm and UTC1 (conversion by `spiceypy`)
- Flags bad fits more or less consistently
- **Lien:**
    - option to fit multiple features (likely best as a class method given the way the code is currently set up--goal is to avoid reading in a given file multiple times)

### Phase Echo Modeling
    Contact: R Maguire
- Stored in `./rss_ringoccs/diffrec/advanced_tools.py` 
- Forward-models diffraction of the uplink calibrating H2 maser to predict location and shape of phase echo in the post-USO profiles
- Requires `Geometry`
- Customizable for particular types of features (edges, wells, hats) and certain radius ranges

## Uranus Pipeline
- Contained in `./rss_uringoccs/`
- **Goal:** reproduce original 1989 processing of Voyager 2 1986 Uranus ring occultation experiments and achieve highest resolution profiles of Uranian rings to date (only possible with a radio occultation experiment).
- Pipeline processing follows MTR86 and Gresh 1989
- Running the pipeline
    - Run from `./rss_uringoccs/pipeline/` with a bash script
    - Parameter file allows user to customize processing before running bash script
- Unique data file reader &mdash; **Contact: J Fong**
    - Reads in unique binary file of 50 kHz Voyager 2 Uranus data digitized by D Simpson. This only contains the real component — imaginary component has been lost!
    - Uses a Hilbert transform to recover imaginary component of the complex signal
- Geometry is computed for each ring frame (feature added by J Fong)
- Completely restructured calibration and DLP steps &mdash; **Contact: S Flury**
    - Raw freq offset is provided, not calculated (50 kHz over a few thousand seconds is computationally taxing)
    - Primary phase steering (as with Saturn) is done for whole data set. Phase detrending function is computed here using a composite Simpson’s rule — it is possible that disagreement with Gresh 1989 results occurs in part because they used a simple Reimann sum to compute the phase detrending function instead of Simpson’s rule. Lien: investigate whether this is in fact the case.
    - DLP resampling occurs at 5 m resolution
    - Secondary phase steering is done by unwrapping phase and “fitting” with a smoothed (convolved with uniform kernel) phase which preserves the diffraction patterns. This was finicky because the phase unwrapping does not handle noise or diffraction patterns well. Lien: investigate possibility of automating this without hardwiring corrections. Perhaps this might make use of the discontinuity-finder to automate the phase unwrapping corrections.
    - Power normalization is performed. Sufficient freespace power exists adjacent to each ring that no checks for gaps is necessary.
- DiffractionCorrection
    - diffrec in rss_ringoccs is generalized to work for any ring occultation and is thus callable from the 
    - A script reads in the Uranus DLP file, builds a DLP object, and calls DiffractionCorrection from rss_ringoccs to do the reconstruction
- **Lien:** 
    - Use Doppler contours to determine location/origin of epsilon ring scattered signal
    - investigate use of `pywavelets` to perform wavelet analysis for Uranian ring system — currently not able to reproduce the results of `wavelet.pro` Morlet wavelet transform due to a missing sine term in the `pywavelets` continuous Morlet wavelet transform

## Other documentation

- readthedocs - https://rss-ringoccs.readthedocs.io/en/latest/
    - Linked to GitHub repository
    - Compiled from sphinx markup in rss_ringoccs doc strings
    - Explains object classes, significant functions and tools, by showing and defining positional and keyword arguments
    - Some math is shown for individual functions as needed
- Quick Start — pdf on GitHub repository
    - Quick outline of all terminal commands for downloading, configuring, and setting up the software
    - Step-by-step outline of using the software without any details of how it works
- User’s Guide — pdf on GitHub repository
    - Conceptual outline
    - Elaboration on installation, configuration, and getting started
    - science/technical justification/explanation for pipeline steps
    - Detailed textual explanation of each step
