from rss_ringoccs import diffrec
from .CSV_tools import ExtractCSVData
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy import interpolate
from .write_output_files import construct_filepath

def cringplots(TauInst, outfile="outfile.pdf"):
    RegDict = [
        ["B1 U1", 74665],
        ["Mimas 4:1", 74900],
        ["B3, B4, U2, R-b", 74949],
        ["B5/U3", 76020],
        ["B6/U7", 76240],
        ["B9/V1", 76430],
        ["76,460", 76460],
        ["B10", 76530],
        ["B11/U8", 76730],
        ["Titan -1:0", 77535],
        ["Titan -1:0 Cont.", 77565],
        ["B13/R-e", 80980],
        ["B14", 81020],
        ["B15/V2", 82000],
        ["B16/R-f", 82050],
        ["B17/R-g", 82200],
        ["B18/R-h", 83625],
        ["B19/R-i", 84635],
        ["B20", 84825],
        ["B21", 84860],
        ["B22", 85110],
        ["B23/R-j", 85450],
        ["B24/U9", 85485],
        ["B25/U4", 85510],
        ["B26/U6", 85530],
        ["B27/R-d", 85680],
        ["B28a", 86415],
        ["B28b", 86445],
        ["B29", 86575],
        ["B31", 86590],
        ["B32/V3", 87185],
        ["Atlas 2:1", 87650],
        ["B35", 88745],
        ["Mimas 6:2", 89890],
        ["Pand 4:2", 89900],
        ["B40", 90280],
    ]

    res = str(TauInst.res)+" km"
    N_Plots = len(RegDict)
    N_Pages = int(N_Plots/8.0)
    
    with PdfPages(outfile) as pdf:
        for i_page in range(N_Pages):
            plt.rc('font', family='serif')
            plt.rc('font', size=10)
            plt.figure(figsize=(8.5, 11))
            plt.suptitle("C-Ring Reconstructions: %s Resolution" % res, size=14)
            gs = gridspec.GridSpec(4, 2, wspace=0.5, hspace=0.5)

            # Plot Normalized Power
            for i_plot in range(8):
                plt.subplot(gs[int(i_plot/2), int(i_plot % 2)])
                plt.tick_params(axis='y', which='both', left=True,
                                right=True, labelleft=True)
                plt.tick_params(axis='x', which='both', bottom=True,
                                top=True, labelbottom=True)
                plt.locator_params(axis='y', nbins=3)
                plt.locator_params(axis='x', nbins=4)
                plt.title(RegDict[int(8*i_page+i_plot)][0])
                plt.ylabel("Normalized Power")
                plt.xlabel("Ring Radius (km)")
                plt.plot(TauInst.rho_km_vals, TauInst.power_vals, 'b')
                rmin = RegDict[int(8*i_page+i_plot)][1] - 15.0
                rmax = RegDict[int(8*i_page+i_plot)][1] + 15.0
                nmin = np.min((TauInst.rho_km_vals >= rmin).nonzero())
                nmax = np.max((TauInst.rho_km_vals <= rmax).nonzero())
                p = TauInst.power_vals[nmin: nmax+1]
                ymin = np.min(p)*0.98
                ymax = np.max(p)*1.02
                plt.xlim(rmin, rmax)
                plt.ylim(ymin, ymax)

            pdf.savefig(bbox_inches="tight", pad_inches=1)
            plt.close()
        
        if ((N_Plots % 8) != 0):
            plt.rc('font', family='serif')
            plt.rc('font', size=10)
            plt.figure(figsize=(8.5, 11))
            plt.suptitle("C-Ring Reconstructions: %s Resolution" % res, size=14)
            gs = gridspec.GridSpec(4, 2, wspace=0.5, hspace=0.5)
            for i_plot in range(int(N_Plots % 8)):
                plt.subplot(gs[int(i_plot/2), int(i_plot % 2)])
                plt.tick_params(axis='y', which='both', left=True,
                                right=True, labelleft=True)
                plt.tick_params(axis='x', which='both', bottom=True,
                                top=True, labelbottom=True)
                plt.locator_params(axis='y', nbins=3)
                plt.locator_params(axis='x', nbins=4)
                plt.title(RegDict[int(8*(N_Pages-1)+i_plot)][0])
                plt.ylabel("Normalized Power")
                plt.xlabel("Ring Radius (km)")
                plt.plot(TauInst.rho_km_vals, TauInst.power_vals, 'b')
                rmin = RegDict[int(8*(N_Pages-1)+i_plot)][1] - 15.0
                rmax = RegDict[int(8*(N_Pages-1)+i_plot)][1] + 15.0
                nmin = np.min((TauInst.rho_km_vals >= rmin).nonzero())
                nmax = np.max((TauInst.rho_km_vals <= rmax).nonzero())
                p = TauInst.power_vals[nmin: nmax+1]
                ymin = np.min(p)*0.98
                ymax = np.max(p)*1.02
                plt.xlim(rmin, rmax)
                plt.ylim(ymin, ymax)

            pdf.savefig(bbox_inches="tight", pad_inches=1)
            plt.close()
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
    nstart = np.max((rec.rho_km_vals == rmax).nonzero())
    nend = np.min((rec.rho_km_vals == rmin).nonzero())

    with PdfPages(outfile) as pdf:
        plt.rc('font', family='serif')
        plt.rc('font', size=10)
        plt.figure(figsize=(8.5, 11))
        plt.suptitle("TAU Parameters", size=14)
        gs = gridspec.GridSpec(3, 2, wspace=0.0, hspace=0.0)

        # Plot Normalized Power
        plt.subplot(gs[0, 0])
        plt.tick_params(axis='y', which='both', left=True,
                        right=False, labelleft=True)
        plt.tick_params(axis='x', which='both', bottom=False,
                        top=True, labelbottom=False)
        plt.title("Comparison Plots")
        plt.ylabel('Normalized Power')
        plt.plot(rec.rho_km_vals, rec.power_vals, 'b')
        plt.plot(data.tau_rho, data.power_vals, 'g')
        plt.xlim(rmin, rmax)
        plt.ylim(-0.1, np.max([np.max(rec.power_vals),
                              np.max(data.power_vals)])*1.1)

        # Plot Difference in Normalized Power
        plt.subplot(gs[0, 1])
        interp = interpolate.interp1d(data.tau_rho, data.power_vals)
        diff = rec.power_vals - interp(rec.rho_km_vals)
        plt.tick_params(axis='y', which='both', left=False,
                        right=True, labelleft=False, labelright=True)
        plt.tick_params(axis='x', which='both', bottom=False,
                        top=True, labelbottom=False)
        plt.title("Difference Plots")
        plt.plot(rec.rho_km_vals, diff, 'r')
        plt.xlim(rmin, rmax)
        plt.ylim(np.min(diff), np.max(diff))

        # Plot Normal Optical Depth
        plt.subplot(gs[1, 0])
        plt.tick_params(axis='y', which='both', left=True, right=False,
                        labelleft=True, labelright=False)
        plt.tick_params(axis='x', which='both', bottom=False, top=False,
                        labelbottom=False, labeltop=False)
        plt.ylabel('Normal Optical Depth')
        plt.plot(rec.rho_km_vals, rec.tau_vals, 'b')
        plt.plot(data.tau_rho, data.tau_vals, 'g')
        plt.xlim(rmin, rmax)
        plt.ylim(-0.1, np.max([np.max(rec.tau_vals),
                              np.max(data.tau_vals)])+0.2)

        # Plot Difference in Normal Optical Depth
        plt.subplot(gs[1, 1])
        plt.xlim(rmin, rmax)
        interp = interpolate.interp1d(data.tau_rho, data.tau_vals)
        diff = rec.tau_vals - interp(rec.rho_km_vals)
        plt.tick_params(axis='y', which='both', left=False,
                        right=True, labelleft=False, labelright=True)
        plt.tick_params(axis='x', which='both', bottom=False,
                        top=False, labelbottom=False)
        plt.plot(rec.rho_km_vals, diff, 'r')
        plt.ylim(np.min(diff), np.max(diff))

        # Plot Residual Phase (Radians)
        ymax = np.max([np.max(rec.phase_vals), np.max(data.phase_vals)])+0.1
        ymin = np.min([np.min(rec.phase_vals), np.min(data.phase_vals)])-0.1
        plt.subplot(gs[2, 0])
        plt.tick_params(axis='y', which='both', left=True, right=False,
                        labelleft=True, labelright=False)
        plt.tick_params(axis='x', which='both', bottom=True, top=False,
                        labelbottom=True, labeltop=False)
        plt.xlabel("Ring Radius (km)")
        plt.ylabel('Residual Phase (Radians)')
        plt.plot(rec.rho_km_vals, rec.phase_vals, 'b', label="rss_ringoccs")
        plt.plot(data.tau_rho, data.phase_vals, 'g', label="PDS")
        plt.legend(loc='upper right', bbox_to_anchor=(1.9, 3.4))
        plt.xlim(rmin, rmax)
        plt.ylim(ymin, ymax)

        # Plot Difference in Residual Phase (Radians)
        plt.subplot(gs[2, 1])
        plt.xlim(rmin, rmax)
        interp = interpolate.interp1d(data.tau_rho, data.phase_vals)
        diff = rec.phase_vals - interp(rec.rho_km_vals)
        plt.tick_params(axis='y', which='both', left=False,
                        right=True, labelleft=False, labelright=True)
        plt.tick_params(axis='x', which='both', bottom=True,
                        top=False, labelbottom=True)
        plt.xlabel("Ring Radius (km)")
        plt.plot(rec.rho_km_vals, diff, 'r')

        pdf.savefig(bbox_inches="tight", pad_inches=1)
        plt.close()
