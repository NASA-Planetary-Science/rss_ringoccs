"""
Purpose: Provide all inputs to batch run, to be imported by e2e_batch.py
"""

### Batch inputs
with16 = False                      # Process 16 kHz files (use False to skip)
mpath = '../data/'                  # Path to location of data files
rsr_file_list = '../tables/rsr_1kHz_files_before_USO_failure.txt'

### Global inputs
verbose    = False                  # Print processing steps to terminal
write_file = True                   # Write output data and label files

### RSRReader
decimate_16khz_to_1khz = False       # Decimate 16 kHz rsr file to 1 kHz

### Geometry
kernels = '../tables/e2e_kernels.ker'  # Path to meta-kernel or list of paths to
                                       #       necessary kernel files
planet = 'Saturn'                      # Name of target planet
spacecraft = 'Cassini'                 # Name of target spacecraft
pt_per_sec = 1.0                       # Number of points per second calculated
                                       #       for each parameter

### Calibration
dt_cal = 1.0                        # Time spacing in seconds between points
pnf_order = 3                       # Power normalization fit order
interact = False                    # Manually update power normalization fit

### DiffractionLimitedProfile
dr_km_desired = 0.25                # Radial spacing between points (or
                                    #       half of desired DLP resolution)
profile_range = [65000., 150000.]   # Lower and upper radial bounds to the 
                                    #       output diffraction profile

### DiffractionCorrection
res_km = 1.0                        # Reconstruction resolution
res_factor = 0.75                   # Factor to be multiplied to resolution
                                    #       to follow MTR86 definition
inversion_range = profile_range     # Lower and upper radial bounds to inversion
psitype = 'Fresnel4'                # Psi type
wtype = 'kbmd20'                    # Window type
sigma = 2.e-13                      # Allen deviation
fwd = False                         # Compute forward model
norm = True                         # Normalize reconstructed complex
                                    #       transmittance by window width
bfac = False                        # Use input sigma in window calculation
