"""
e2e_run.py

Purpose:
    Example rss_ringoccs package script to produce a 1-km resolution
    profile of the Huygens ringlet for Rev007 E X43.

Notes:
    #. This program assumes that all relevant files (rsr_file and
        kernels) have already been downloaded. See p6-7 of
        rss_ringoccs: User's Guide for details.
"""
import e2e_run_args as args
import sys
import time
import matplotlib.pyplot as plt
import numpy as np

sys.path.append('../')
import rss_ringoccs as rss
sys.path.remove('../')

start_time = time.time()

# Create instance with rsr file contents
rsr_inst = rss.rsr_reader.RSRReader(args.rsr_file, verbose=args.verbose,
        decimate_16khz_to_1khz=args.decimate_16khz_to_1khz)

# Create instance with geometry parameterst
geo_inst = rss.occgeo.Geometry(rsr_inst, args.planet, args.spacecraft,
        args.kernels, verbose=args.verbose, write_file=args.write_file)

# Create instance with calibrated data
cal_inst = rss.calibration.Calibration(rsr_inst, geo_inst, verbose=args.verbose,
        write_file=args.write_file, fof_order=args.fof_order, 
        pnf_order=args.pnf_order, pnf_fittype=args.pnf_fittype,
        interact=args.interact)

# Create instance with diffraction-limited profile and other
#   inputs needed for diffraction correction
dlp_inst_ing, dlp_inst_egr = (
        rss.calibration.DiffractionLimitedProfile.create_dlps(
            rsr_inst, geo_inst, cal_inst, args.dr_km_desired, 
            profile_range = args.profile_range,
            write_file=args.write_file, verbose=args.verbose))

# Invert profile for full occultation
if dlp_inst_ing is not None:
    tau_inst = rss.diffrec.DiffractionCorrection(dlp_inst_ing, args.res_km,
            rng=args.inversion_range, res_factor=args.res_factor,
            psitype=args.psitype, wtype=args.wtype, fwd=args.fwd,
            norm=args.norm, bfac=args.bfac, write_file=args.write_file,
            verbose=args.verbose)
if dlp_inst_egr is not None:
    tau_inst = rss.diffrec.DiffractionCorrection(dlp_inst_egr, args.res_km,
            rng=args.inversion_range, res_factor=args.res_factor,
            psitype=args.psitype, wtype=args.wtype, fwd=args.fwd,
            norm=args.norm, bfac=args.bfac, write_file=args.write_file,
            verbose=args.verbose)


# Plot Huygens ringlet
#   first row of three: uncorrected optical depth, power, and phase
#   second row of three: corrected optical depth, power, and phase
fig, axes = plt.subplots(2,3, figsize=(8,5))

rho_km = tau_inst.rho_km_vals - args.feature_km
recon_power = tau_inst.power_vals
raw_power = tau_inst.p_norm_vals

raw_tau = -tau_inst.mu_vals * np.log(tau_inst.p_norm_vals)
recon_tau = tau_inst.tau_vals

raw_phase = tau_inst.phase_rad_vals
recon_phase = tau_inst.phase_vals
xtitle = '$\\rho$ - ' + str(args.feature_km) + ' (km)'

fig.suptitle('Rev007E X43 ' + args.feature_name + ' at ' + str(args.res_km) + ' km resolution', fontsize=14)

axes[0,0].plot(rho_km, raw_power)
axes[0,0].set_ylabel('Raw Power')
axes[0,0].set_xlabel(xtitle)

axes[0,1].plot(rho_km, raw_tau)
axes[0,1].set_ylabel('Raw Optical Depth')
axes[0,1].set_xlabel(xtitle)

axes[0,2].plot(rho_km, raw_phase)
axes[0,2].set_ylabel('Raw Phase (rad)')
axes[0,2].set_xlabel(xtitle)

axes[1,0].plot(rho_km, recon_power)
axes[1,0].set_ylabel('Reconstructed Power')
axes[1,0].set_xlabel(xtitle)

axes[1,1].plot(rho_km, recon_tau)
axes[1,1].set_ylabel('Reconstructed Optical Depth')
axes[1,1].set_xlabel(xtitle)

axes[1,2].plot(rho_km, recon_phase)
axes[1,2].set_ylabel('Reconstructed Phase (rad)')
axes[1,2].set_xlabel(xtitle)

end_time = time.time()
plt.tight_layout()
plt.subplots_adjust(top=0.93)

print('Total run time (min): ', str((end_time-start_time)/60.))

plt.show()

"""
Revisions:
"""
