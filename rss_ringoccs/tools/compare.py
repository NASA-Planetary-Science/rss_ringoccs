from rss_ringoccs import diffrec
from .CSV_tools import ExtractCSVData
import spiceypy as spice
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy import interpolate
from .write_output_files import construct_filepath


# title_out, outdir_out = construct_filepath(rev_info, filetyp)

def compare(NormDiff, geo, cal, dlp, tau, outfile, res=0.75, rng="Maxwell",
            wtype="kbmd20", norm=True, bfac=True, sigma=2.e-13, verbose=True,
            psitype="Fresnel8"):
    data = ExtractCSVData(geo, cal, dlp, tau=tau, verbose=verbose)
    rec = diffrec.DiffractionCorrection(NormDiff, res, rng=rng, wtype=wtype,
                                        fwd=True, norm=norm, bfac=bfac,
                                        sigma=sigma, verbose=verbose,
                                        psitype=psitype, write_file=False)
    rmin = np.min(rec.rho_km_vals)
    rmax = np.max(rec.rho_km_vals)

    with PdfPages(outfile) as pdf:
        plt.rc('font', family='serif')
        plt.rc('font', size=10)
        plt.figure(figsize=(8.5, 11))
        plt.suptitle("TAU Parameters", size=14)
        gs = gridspec.GridSpec(3, 2, wspace=0.0, hspace=0.0)

        # Plot Normalized Power
        plt.subplot(gs[0, 0])
        plt.xlim(rmin, rmax)
        plt.tick_params(axis='y', which='both', left=True,
                        right=True, labelleft=True)
        plt.tick_params(axis='x', which='both', bottom=False,
                        top=True, labelbottom=False)
        plt.ylabel('Normalized Power')
        plt.plot(rec.rho_km_vals, rec.power_vals, label="rss_ringoccs")
        plt.plot(data.tau_rho, data.power_vals, label="PDS")
        plt.legend()

        # Plot Difference in Normalized Power
        plt.subplot(gs[0, 1])
        plt.xlim(rmin, rmax)
        interp = interpolate.interp1d(data.tau_rho, data.power_vals)
        power = interp(rec.rho_km_vals)
        plt.tick_params(axis='y', which='both', left=True,
                        right=True, labelleft=True)
        plt.tick_params(axis='x', which='both', bottom=False,
                        top=True, labelbottom=False)
        plt.ylabel('Difference in Power')
        plt.plot(rec.rho_km_vals, rec.power_vals-power)
        plt.legend()

        plt.subplot(gs[1,0])
        plt.xlim(rmin, rmax)
        plt.tick_params(axis='y', which='both', left=True, right=True,
                        labelleft=True, labelright=False)
        plt.tick_params(axis='x', which='both', bottom=False, top=False,
                        labelbottom=False, labeltop=False)
        plt.ylabel('Normal Optical Depth')
        plt.plot(rec.rho_km_vals, rec.tau_vals, label="rss_ringoccs")
        plt.plot(data.tau_rho, data.tau_vals, label="PDS")
        plt.legend()
        plt.subplot(gs[2,0])
        plt.xlim(rmin, rmax)
        plt.tick_params(axis='y', which='both', left=True, right=True,
                        labelleft=True, labelright=False)
        plt.tick_params(axis='x', which='both', bottom=True, top=False,
                        labelbottom=True, labeltop=False)
        plt.xlabel("Ring Radius (km)")
        plt.ylabel('Residual Phase (Radians)')
        plt.plot(rec.rho_km_vals, rec.phase_vals, label="rss_ringoccs")
        plt.plot(data.tau_rho, data.phase_vals, label="PDS")
        plt.legend()

        pdf.savefig(bbox_inches="tight",pad_inches=1)
        plt.close()


# B_rad_vals
# D_km_vals
# F_km_val
# f_sky_hz_vals
# p_norm_fwd_vals
# p_norm_vals
# phase_fwd_vals
# phase_rad_vals
# phase_vals
# phi_rad_vals
# phi_rl_rad_vals
# power_vals
# raw_tau_threshold_vals
# rho_dot_kms_vals
# rho_km_vals
# tau_threshold_vals
# tau_vals
