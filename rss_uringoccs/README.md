# rss_uringoccs
rss_uringoccs is a suite of open-source Python-based analysis tools for the Voyager 2 Uranus Radio Science (RSS) ring occultations. It was developed by Team Cassini at Wellesley College (Sophia Flury, Jolene Fong, Ryan Maguire, and Glenn Steranka) under the direction of Richard French, Cassini RSS Team Leader. rss_uringoccs is based on rss_ringoccs, which was developed to process RSS Saturn ring observations from Cassini and Voyager 2. 

## Introduction
In January of 1986, Voyager 2 performed a radio occultation experiment of the Uranus ring system. The purpose of rss_uringoccs is to enable scientists to produce "on demand" radial optical depth profiles of the Uranian rings from the raw RSS data, without requiring a deep familiarity with the complex processing steps involved in calibrating the data and correcting for the effects of diffraction. The code and algorithms are documented, providing a starting point for users who wish to test, refine, or optimize the straightforward methods we have employed. Our emphasis has been on clarity, sometimes at the expense of programming efficiency and execution time. rss_uringoccs does an excellent job of reproducing existing 50 m-resolution RSS processed ring occultation data available on NASA's Planetary Data System (PDS) at the Ring-Moon Systems Node, but we make no claim to having achieved the state-of-the-art in every respect. We encourage users to augment our algorithms and to report on those improvements, so that they can be  incorporated in future editions of rss_uringoccs.

The original processing of the Voyager 2 Uranus RSS ring occultation observations was performed by Gresh et al. (Icarus 78, 131-168, 1989) [1989Icar...78..131G](https://www.sciencedirect.com/science/article/abs/pii/0019103589900742?via%3Dihub). Their methodology follows from the diffraction reconstruction technique presented by Marouf, Tyler, and Rosen (Icarus 68, 120-166, 1986) [1986Icar...68..120M](https://www.sciencedirect.com/science/article/abs/pii/0019103586900783?via%3Dihub).

## Installation and Documentation
The fundamental instructions for using rss_uringoccs are discussed below. For details regarding dependencies, architecture, and methods, we refer the user to the rss_ringoccs documentation on GitHub and readthedocs.io.

Detailed installation instructions and full documentation for rss_ringoccs are contained in https://github.com/NASA-Planetary-Science/rss_ringoccs/tree/master/docs/rss_ringoccs_User_Guide.pdf.
The dependencies are identical, and the directory structure and file naming conventions are very similar.

Release notes for rss_ringoccs are contained in https://github.com/NASA-Planetary-Science/rss_ringoccs/blob/master/ReleaseNotes.md

Source code documentation for rss_ringoccs is found at https://rss-ringoccs.readthedocs.io/en/master/

## Directory Structure and Contents
Below is a brief outline of the top-level directory structure and contents of the rss_uringoccs package.

|DIRECTORY								| CONTAINS																																															|
|-------------------------|-------------------------------------------------------------------------------------------------------|
|rss_uringoccs/data/			| this is where the raw data files should be stored	for the pipeline to access and process							|
|rss_uringoccs/examples/	| this contains example scripts for reading in and plotting the output science products									|
|rss_uringoccs/kernels/		| this contains the reference NAIF/SPICE kernels for computing the occultation geometry									|
|rss_uringoccs/output/		| this contains the subdirectories of output science products as processed by the pipeline							|
|rss_uringoccs/pipeline/	| the pipeline scripts which process the Voyager 2 Uranus data with an adapted version of rss_ringoccs	|
|rss_uringoccs/tables/		| the reference files for processing the raw data, including ancillary frequency offsets and kernels		|

## Data and reference files
The requisite raw data and reference kernel files have not yet been archived to NASA's PDS, but they are available upon request from Richard French (rfrench@wellesley.edu). The retrieved compressed file should be unzipped in the `rss_uringoccs` directory. This will create a new directory within `rss_uringoccs` called `uring_files`. Within `uring_files` there will be two directories: `kernels` and `data`. The `kernels` directory containa NAIF/SPICE kernel files necessary for computing the geometry of the Uranus ring system, a prerequisite for reconstructing the ring profiles. The `data` directory contains binary files containing raw Voyager 2 Uranus rings occultation data recorded at 50 kHz.

Additional ancillary files are included in the `rss_uringoccs/tables/` directory, including the offset and sky frequencies computed from the raw 50 kHz data.

## Running rss_uringoccs
Unlike the rss_ringoccs software from which this was developed, the rss_uringoccs pipeline requires running multiple scripts. This difference in organization is in part due to ring reference frames in the occultation geometry, high-resolution raw data, a large frequency offset in the tens of kHz, and an additional phase steering step. The goal of distributing the pipeline over multiple scripts is to ultimately reduce the computation time for the user. Each script is designed to process all 9 ring profiles for both ingress and egress and thus needs to be run only once for a given set of inputs.

Input values and options are managed in the `rss_uringoccs/pipeline/pipeline_params.py` script. This is the only file users might wish  to edit. It contains the arguments for input files (e.g., kernel files, GEO and CAL files, etc), DLP sampling, and reconstruction options (e.g., window, processing resolution, etc.).

To run the entire pipeline, navigate to the `rss_uringoccs/pipeline/` directory and run the `uring_pipeline.sh` shell script with the appropriate permissions:
```
cd rss_uringoccs/pipeline/
chmod +x uring_pipeline.sh
bash uring_pipeline.sh
```
This script will call each of the pipeline scripts --`geocal.py`, `dlp.py`, and `reconstruction.py`-- in succession.

This will produce reduced science data files in the `rss_uringoccs/output/` directory under the appropriate profile direction and ring subdirectories. For instance, to find the science products for the delta ring egress profile, one would look in `rss_uringoccs/output/E/D/`. File naming convention is
```
VGR2_X43_P_URING_R_TAU_XXXXXM_YYYYMMDD_SSSS.TAB
```
where `P` is the profile direction, `R` is the ring identifier, `XXXXX` is the resolution in meters, `YYYYMMDD` is the processing date, and `SSSS` is the serial number for processing results on that date. For example, for the first set of delta ring egress results processed on May 31, 2019 at 20 m resolution, the reconstructed profile would be stored in
```
rss_uringoccs/output/E/D/VGR2_X43_E_URING_D_TAU_00020M_20190531_0001.TAB
```
with a reference documentation file
```
rss_uringoccs/output/E/D/VGR2_X43_E_URING_D_TAU_00020M_20190531_0001.LBL
```

Also included in the `rss_uringoccs/pipeline/` directory is a script `scatter.py` for processing the incoherent signal accompanying the carrier signal from Voyager 2. This script was used to reproduce the epsilon ring scattered signal results in Gresh et al. (1989). To read and/or visualize scattered signal results, see the `rss_uringoccs/examples/spectrogram_plot_example.py` script. The scripts `profile_plots.py` and `gallery_plots.py` in `rss_uringoccs/examples/` are also provided as a guide for visualizing the reconstructed optical depth profiles for each Uranian ring. For the purpose of validation, we also include `compare_gresh.py` in this directory to compare a user's reconstructed optical depth profiles to those from Gresh et al. (1989).

We have benchmarked the `uring_pipeline.sh` script at 30 min run time to process all 9 rings in both ingress and egress at 50 m resolution reconstructed profiles. However, this benchmark time will likely vary from system to system and could exceed more than an hour of computation time.

## How to (Get) Help
If you have trouble with installation or execution of the rss_uringoccs package, we encourage you to post a issue to https://github.com/NASA-Planetary-Science/rss_ringoccs/issues. We will attempt to respond promptly, and ther users will benefit. Alternatively, you can write email directly to Richard French: rfrench_at_wellesley.edu.
## Citing rss_uringoccs
If you use rss_uringoccs as the basis of a publication, please consider
citing rss_ringoccs using the [DOI:10.5281/zenodo.2548947](https://doi.org/10.5281/zenodo.2557755)

## Acknowledgements
Development of the Uranus radio occultation processing scripts in rss_uringoccs was supported by NASA program NNH14ZDA001N-PDART; rss_ringoccs was developed with finding provided by the Cassini project as part of its mission closeout efforts. We are especially grateful to Donna Gresh and her colleagues for the initial processing of the Voyager 2 Uranus observations, to Dick Simpson (NASA/PDS) for tracking down the elusive Voyager 2 RSS ring occultation observations, and to Cassini RSS Team Member Essam Marouf for developing the diffraction reconstruction technique that underlies this work.
