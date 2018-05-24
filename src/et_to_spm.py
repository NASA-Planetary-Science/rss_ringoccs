
"""

et_to_spm.py

Purpose: Convert ephemeris time to seconds past midnight (SPM).

NOTE (GJS): Input ephemeris time must be increasing

Revisions:
	2018 Mar 13 - jfong - original
	2018 May 23 - gsteranka - ibrk[0] --> int(ibrk[0])
"""

import numpy as np
import spiceypy as spice

def et_to_spm(et_vals, kernels=None):
	"""Function to convert ephemeris time to SPM"""

	npts = len(et_vals)

	spm_vals = []

	for n in range(npts):
		et = et_vals[n]

		# Convert ET to UTC string
		utc_str = spice.et2utc(et, "ISOD", 16)
		
		# Extract hour, minute, and second from UTC string
		hour = float(utc_str[9:11])
		minute = float(utc_str[12:14])
		second = float(utc_str[15:])
		spm = hour*3600. + minute*60. + second
		spm_vals.append(spm)

	# Check if event goes through midnight
	dr_vals = spm_vals - np.roll(spm_vals, 1)
	ibrk = np.argwhere(dr_vals[1:] < 0)

	if ibrk.size > 0:
		# CHANGE (GJS): Need to "int" this so it's not an array
		ind = int(ibrk[0]) + 1
		spm_brk = spm_vals[ind]
		spm_vals_cont = np.zeros(npts)
		spm_vals_cont[0:ind] = spm_vals[0:ind]
		spm_vals_cont[ind+1:] = spm_vals[ind+1:] + spm_brk
		spm_vals = spm_vals_cont
	
	return spm_vals


