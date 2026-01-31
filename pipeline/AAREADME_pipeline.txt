rss_ringoccs/pipeline/AAREADME_pipeline.txt

rss_ringoccs/pipeline/AAREADME_pipeline.txt 

Revisions:
  2026 Jan 30 - rfrench 

The rss_ringoccs/pipeline/ directory contains automated scripts 
to perform end-to-end processing of raw Cassini RSS ring occultation 
data to produce diffraction-reconstructed radial optical depth profiles.

The initial steps are to obtain the following data files from the web:

- kernel files that include planetary and spacecraft ephemerides, 
  planetary constants and DSN coordinates

  The following script retrieves the necessary kernels from NAIF website and 
  populates the rss_ringoccs/kernels/ directory

  get_kernels.sh  

  The directory structure is the same as on the NAIF website:

	rss_ringoccs/kernels
	├── local
	├── naif
	│   ├── CASSINI
	│   │   └── kernels
	│   │       ├── fk
	│   │       │   ├── prelim
	│   │       │   └── stations
	│   │       │       └── prelim
	│   │       ├── lsk
	│   │       ├── pck
	│   │       └── spk
	│   └── generic_kernels
	│       ├── fk
	│       │   └── stations
	│       │       ├── a_old_versions
	│       │       └── prelim
	│       └── spk
	│           ├── planets
	│           └── stations


- raw (RSR) data files at 1 kHz and optionally at 16 kHz for 
  higher-resolution diffraction reconstructions

 The following scripts retrieve the raw RSR data files. 

 get_1kHz_rsr_files_preUSOfailure.sh 

 get_16kHz_rsr_files_preUSOfailure.sh

 The following script retrieves 1 kHz RSR files for observations after Rev133.
 They are not officially supported by rss_ringoccs because of phase instability.
 However, they are useful for diffraction-limited studies and selected isolated
 features that can be successfully diffraction-reconstructed.

 get_1kHz_rsr_files_postUSOfailure.sh

 The retrieved data files have the following directory structure to match
 the source directory structure on the PDS website:

	rss_ringoccs/data
	├── co-s-rss-1-sroc1-v10
	│   ├── cors_0105
	│   │   └── sroc1_123
	│   │       └── rsr
	│   ├── cors_0106
	│   │   └── sroc1_123
	│   │       └── rsr
	│   ├── cors_0107
	│   │   └── sroc1_123
	│   │       └── rsr
	│   ├── cors_0108
	│   │   └── sroc1_123
	│   │       └── rsr
	│   ├── cors_0109
	│   │   └── sroc1_141
	│   │       └── rsr


- the CORSS_8001 NASA/PDS/Ring-moon-systems archive of 1 km resolution Cassini 
  RSS ring profiles is used by rss_ringoccs to provide comparisons with
  independently produced radial profiles 

 This archive can be obtained from:

 https://pds-rings.seti.org/holdings/archives-volumes/CORSS_8xxx/CORSS_8001.tar.gz

 The following python program downloads and unpacks this compressed volume if it is not
 already present:

 python3 get_CORSS_8001.py
 
 The directory structure of the extracted archive is shown below:

	rss_ringoccs/data/CORSS_8001
	├── browse
	├── catalog
	├── data
	│   ├── Rev007
	│   │   ├── Rev007E
	│   │   │   ├── Rev007E_RSS_2005_123_K34_E
	│   │   │   ├── Rev007E_RSS_2005_123_S43_E
	│   │   │   ├── Rev007E_RSS_2005_123_X34_E
	│   │   │   └── Rev007E_RSS_2005_123_X43_E
 

Once all of these files are in place, the user can execute the following command
to produce 1km resolution optical depth profiles (TAU files), along with the
auxiliary GEO, CAL, and DLP files, for the full RSS data set from Rev007 to Rev137:


The script:
	e2e_batch_1km_20260130.sh
executes the following command:
	echo yes | e2e_batch_1km_20260130.py

'echo yes' is required to provide an automated response to 
queries to the user by rss_ringoccs

The run time for e2e_batch_1km_20260130.sh is 

when run on the following system:

Mac Studio 2022 Apple M1 Ultra MacOS Tahoe 26.2 128 GB RAM python3.14


The previous procedure uses only a single core. Execution time can be 
greatly reduced by using multiprocessing in Python. 
A Mac Studio 2022 has 24 logical cores, all of which can be used
simultaneously. The main restriction is that rss_ringoccs can use
quite a bit of RAM for high-resolution reconstructions across a wide
ring region.


The following Python procedure uses multiprocessing to read 1 kHz
raw RSR files and produce GEO, CAL, and DLP (diffraction-limited profiles)
at a variety of resolutions for all observations prior to the USO failur.

python3 batch_parallel_raw_to_DLP_1kHz_to_50m_20260130.py

The run time for python3 batch_parallel_raw_to_DLP_1kHz_to_50m_20260130.py
is xx hours when run on the following system:

Mac Studio 2022 Apple M1 Ultra MacOS Tahoe 26.2 128 GB RAM python3.14

Once the GEO, CAL, and DLP files have been produced, they can be used for
a wide variety of diffraction construction request. They need to be 
produced only once.

The following procedures generate TAU files from the GEO, CAL, and DLP files.

First, the following Jupyter notebook (or, alernatively, the corresponding *.py
file) specifies the specific event files to be processed, and the desired
radial range and resolution of the diffraction reconstruced TAU files.

construct_batch_parallel_CSVtoTAU_20260130.ipynb
construct_batch_parallel_CSVtoTAU_20260130.py

The method used is to populate the following template file:

batch_parallel_CSVtoTAU_template-20260130a.py

and to replace $-delimited keywords in the template with specific values
specified in construct_batch_parallel_CSVtoTAU_20260130.

Each newly-created python procedure is saved to 

rss_ringoccs/pipeline/runs_20260130a/

The output of batch_parallel_CSVtoTAU_template-20260130a.py lists the conditions
for each job. For example:

runs_20260130a/batch_parallel_CSVtoTAU-20260130a-0000.py 1.0 Rev007 I X26 NewtonFilon12 C + inner B ring
runs_20260130a/batch_parallel_CSVtoTAU-20260130a-0001.py 1.0 Rev007 I X43 NewtonFilon12 C + inner B ring
runs_20260130a/batch_parallel_CSVtoTAU-20260130a-0002.py 1.0 Rev007 I K26 NewtonFilon12 C + inner B ring
runs_20260130a/batch_parallel_CSVtoTAU-20260130a-0003.py 1.0 Rev007 I S43 NewtonFilon12 C + inner B ring

the '1.0' means 1 km reesolution
the event is Rev007_X26I (X band, DSS-26)
the quadrature method (psitype) is Newton Filon 12
The ring region is C + inner B ring (74400-97000 km)

The list of parallel-processing jobs to be run is contained in 

batch_parallel_CSVtoTAU-20260130a.sh

The contents look like this:

python3 runs_20260130a/batch_parallel_CSVtoTAU-20260130a-0000.py
python3 runs_20260130a/batch_parallel_CSVtoTAU-20260130a-0001.py
python3 runs_20260130a/batch_parallel_CSVtoTAU-20260130a-0002.py
python3 runs_20260130a/batch_parallel_CSVtoTAU-20260130a-0003.py
python3 runs_20260130a/batch_parallel_CSVtoTAU-20260130a-0004.py

To run all the jobs in succession, type

batch_parallel_CSVtoTAU-20260130a.sh

To run a specific job, type the specific entry you'd like to run.
For example:

python3 runs_20260130a/batch_parallel_CSVtoTAU-20260130a-0000.py
