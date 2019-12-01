import numpy
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
sys.path.append("../")
from rss_ringoccs.tools.CSV_tools import ExtractCSVData
import rss_ringoccs.diffrec.advanced_tools as at
import rss_ringoccs.diffrec.diffraction_correction as dc


def test_CompareTAU():

    # Location of the data for Rev133E.
    geo = "./Test_Data/Rev133E_X43_CRing_GEO.TAB"
    cal = "./Test_Data/Rev133E_X43_CRing_CAL.TAB"
    dlp = "./Test_Data/Rev133E_X43_CRing_DLP_500M.TAB"
    tau = "./Test_Data/Rev133E_X43_CRing_TAU_1000M.TAB"

    # Name of the rev.
    rev = "Rev133E X43"

    # Resolution for reconstruction (in kilometers).
    res = 1.0

    # Resolution factor for reconstruction.
    res_f = 0.75

    # Range of the reconstruction (the location of the C ring).
    rng = [74500, 90500]

    # Select window type for reconstruction.
    wtype="kbmd20"

    # The approximation to the Fresnel kernel that is used.
    psitype="full"

    # Perturb the Fresnel kernel by a small cubic polynomial.
    perturb = [0,0,0,0.3,0]

    # Name of the output PDF file.
    outfile = "./Figures/CompareTAU_Test.pdf"

    # Run diffraction correction on the input data.
    TauInst = at.CompareTau(geo, cal, dlp, tau, res, rng=rng, wtype=wtype,
                            verbose=False, psitype=psitype, res_factor=res_f,
                            perturb=perturb)

    # List of regions to plot within the C Ring.
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
        ["C Ring Ripples", 77725],
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

    # Convert to a string for plots.
    res = "%s km" % str(TauInst.res)

    # Total number of plots, and total number of pages.
    N_Plots = len(RegDict)
    N_Pages = int(N_Plots/8.0)

    with PdfPages(outfile) as pdf:
        for i_page in range(N_Pages):
            # Font and configuration for plots.
            plt.rc('font', family='serif')
            plt.rc('font', size=10)
            plt.figure(figsize=(8.5, 11))
            plt.suptitle("%s C-Ring Reconstructions: %s Resolution"
                         % (rev, res), size=14)
            
            # Use gridspec to create 4-by-2 plots.
            gs = gridspec.GridSpec(4, 2, wspace=0.5, hspace=0.5)

            # Plot Normalized Power
            for i_plot in range(8):
                plt.subplot(gs[int(i_plot/2), int(i_plot % 2)])

                # Setup parameters for x and y axes.
                plt.tick_params(axis='y', which='both', left=True,
                                right=True, labelleft=True)
                plt.tick_params(axis='x', which='both', bottom=True,
                                top=True, labelbottom=True)

                # Number of ticks on the x and y axes.
                plt.locator_params(axis='y', nbins=3)
                plt.locator_params(axis='x', nbins=4)

                # Set the title, and labels for both axes.
                plt.title(RegDict[int(8*i_page+i_plot)][0])
                plt.ylabel("Normalized Optical Depth")
                plt.xlabel("Ring Radius (km)")

                # Plot the plots.
                plt.plot(TauInst.rho_km_vals, TauInst.tau_vals, 'b')
                plt.plot(TauInst.rho_km_vals, TauInst.tau_tau, 'r')

                # Compute the range for the x-axis.
                rmin = RegDict[int(8*i_page+i_plot)][1] - 25.0
                rmax = RegDict[int(8*i_page+i_plot)][1] + 25.0

                # Get the indices that are being plotted.
                nmin = numpy.min((TauInst.rho_km_vals >= rmin).nonzero())
                nmax = numpy.max((TauInst.rho_km_vals <= rmax).nonzero())

                # Add a buffer in the y-axis.
                p = TauInst.tau_vals[nmin: nmax+1]
                ymin = numpy.min(p)*0.98
                ymax = numpy.max(p)*1.02

                # Set the range for the x and y axes.
                plt.xlim(rmin, rmax)
                plt.ylim(ymin, ymax)

            # Add the plot to the output PDF and move on to the next one.
            pdf.savefig(bbox_inches="tight", pad_inches=1)
            plt.close()

        # Special case for if the number of plots is not divisble by eight.
        if ((N_Plots % 8) != 0):
            # Configuration for the plots.
            plt.rc('font', family='serif')
            plt.rc('font', size=10)
            plt.figure(figsize=(8.5, 11))
            plt.suptitle("%s C-Ring Reconstructions: %s Resolution"
                         % (rev, res), size=14)

            # Use gridspec to create 4-by-2 plots.
            gs = gridspec.GridSpec(4, 2, wspace=0.5, hspace=0.5)

            for i_plot in range(int(N_Plots % 8)):
                plt.subplot(gs[int(i_plot/2), int(i_plot % 2)])

                # Tick parameters for x and y axes.
                plt.tick_params(axis='y', which='both', left=True,
                                right=True, labelleft=True)
                plt.tick_params(axis='x', which='both', bottom=True,
                                top=True, labelbottom=True)
                
                # Number of ticks for the x and y axes.
                plt.locator_params(axis='y', nbins=3)
                plt.locator_params(axis='x', nbins=4)

                # Title and labels for axes.
                plt.title(RegDict[int(8*(N_Pages-1)+i_plot)][0])
                plt.ylabel("Normalized Power")
                plt.xlabel("Ring Radius (km)")

                # Plot the plots.
                plt.plot(TauInst.rho_km_vals, TauInst.tau_vals, 'b')
                plt.plot(TauInst.rho_km_vals, TauInst.tau_tau, 'r')

                # Compute the range of the x and y axes.
                rmin = RegDict[int(8*(N_Pages-1)+i_plot)][1] - 15.0
                rmax = RegDict[int(8*(N_Pages-1)+i_plot)][1] + 15.0

                # Get the indices that are being plotted.
                nmin = numpy.min((TauInst.rho_km_vals >= rmin).nonzero())
                nmax = numpy.max((TauInst.rho_km_vals <= rmax).nonzero())

                # Add a buffer to the y-axis.
                p = TauInst.tau_vals[nmin: nmax+1]
                ymin = numpy.min(p)*0.98
                ymax = numpy.max(p)*1.02

                # Set the range for the x and y axes.
                plt.xlim(rmin, rmax)
                plt.ylim(ymin, ymax)

            # Add the final page to the output PDF and close.
            pdf.savefig(bbox_inches="tight", pad_inches=1)
            plt.close()

def test_ModelFromGEO():

    # Location of the data for Rev007E.
    geo = "./Test_Data/Rev007E_X43_Maxwell_GEO.TAB"
    cal = "./Test_Data/Rev007E_X43_Maxwell_CAL.TAB"
    dlp = "./Test_Data/Rev007E_X43_Maxwell_DLP_500M.TAB"

    # Resolution of reconstruction.
    res = 0.2

    # Wavelength of light (in kilometers).
    wav = 3.6e-5

    # Center of the model being analyzed.
    rho = 87500

    # Range of the input data being processed.
    rng = [86000, 89000]

    # Sample size between data points, in kilometers.
    dxs = 0.05

    # Name of the output PDF file.
    outfile = "./Figures/ModelFromGEO_Test.pdf"

    # Some labels for the plots.
    lbl = "Forward Model"
    lbf = "Fresnel Model"
    lba = "Actual Model"
    rlb = "Reconstruction from Forward Model"
    rfm = "Reconstruction from Fresnel Model"

    # Range in the x-axis for the plots.
    rmin = 87400
    rmax = 87600

    # The width for ringlet, gap, and square wave models, in kilometers.
    width = 100

    # Number of waves calculus per point in the square wave model.
    N_Waves = 20

    with PdfPages(outfile) as pdf:

        # Loop over all models.
        for model in ["Ringlet", "Gap", "Right Straight Edge",
                      "Left Straight Edge", "Square Wave", "Delta Impulse"]:

            # Loop over two different opacities.
            for opacity in {0.0, 0.3}:

                # There is no opacity option for square wave and dirac-delta.
                if (opacity != 0) and ((model == "Square Wave") or
                                       (model == "Delta Impulse")):
                    continue

                # Create the modeled data without the Fresnel approximation.
                gdata = at.ModelFromGEO(geo, wav, res, rho, rng=rng,
                                        dx_km_desired=dxs, model=model,
                                        opacity=opacity, verbose=False,
                                        width=width, N_Waves=N_Waves)

                # Create the modeled data with the Fresnel approximation.
                gdataf = at.ModelFromGEO(geo, wav, res, rho, rng=rng,
                                         dx_km_desired=dxs, use_fresnel=True,
                                         model=model, opacity=opacity,
                                         verbose=False, width=width,
                                         N_Waves=N_Waves)

                # Run diffraction correction on both sets of modeled data.
                rec = dc.DiffractionCorrection(gdata, res, rng=rng)
                rec_f = dc.DiffractionCorrection(gdataf, res, rng=rng)

                # Configurations for the plots.
                plt.rc('font', family='serif')
                plt.rc('font', size=10)
                plt.figure(figsize=(8.5, 11))

                # Title for the page.
                plt.suptitle("ModelFromGEO - %s - 200m Model - Opacity: %.1f"
                             % (model, opacity), size=14)

                # Use gridspec to create a 2-by-1 group of plots.
                gs = gridspec.GridSpec(2, 1, hspace=0.0)

                # Select the first plot region.
                plt.subplot(gs[0, 0])

                # Set up the tick parameters for the y axis.
                plt.tick_params(axis='y', which='both', left=True,
                                right=True, labelleft=True)
                plt.locator_params(axis='y', nbins=4)

                # Title for the plot, and tick parameters for the x-axis.
                plt.title("Maxwell Ringlet - Rev007E Geometry Data")
                plt.tick_params(axis='x', which='both', bottom=False,
                                top=True, labelbottom=False)
                plt.locator_params(axis='x', nbins=8)

                # Set the label for the y-axis.
                plt.ylabel("Normalized Power")

                # Plot the modeled power.
                plt.plot(gdataf.rho_km_vals, gdataf.p_norm_vals, 'r', label=lbf)
                plt.plot(gdata.rho_km_vals, gdata.p_norm_vals, 'b', label=lbl)

                # Set x axis range, and put the legend on top of plots.
                plt.xlim(rmin, rmax)
                plt.legend()
    
                if model == "Delta Impulse":
                    plt.ylim(0, 1.2)

                # Select the second plot region and repeat the above steps.
                plt.subplot(gs[1, 0])

                # Set up the tick parameters for the y axis.
                plt.tick_params(axis='y', which='both', left=True,
                                right=True, labelleft=True)
                plt.locator_params(axis='y', nbins=4)

                plt.xlabel("Ring Radius (km)")
                plt.ylabel("Normalized Power")
                plt.tick_params(axis='x', which='both', bottom=True,
                                top=False, labelbottom=True)
                plt.locator_params(axis='x', nbins=8)

                # Plot the reconstructed power.
                plt.plot(gdata.rho_km_vals, gdata.p_norm_actual_vals,
                         'g', label=lba)
                plt.plot(rec_f.rho_km_vals, rec_f.power_vals, 'r', label=rfm)
                plt.plot(rec.rho_km_vals, rec.power_vals, 'b', label=rlb)

                plt.xlim(rmin, rmax)
                plt.legend()

                # Append the plots to the pdf, close, and move on to the next one.
                pdf.savefig(bbox_inches="tight", pad_inches=1)
                plt.close()

        # Test the arbitrary data modeling tool.
        data = ExtractCSVData(geo, cal, dlp, verbose=False)
        rec = dc.DiffractionCorrection(data, 1.0, rng=[85000, 89000])
        gdata = at.ModelFromGEO(geo, wav, res, rho, rng=rng, dx_km_desired=dxs,
                                model="fromdata", data_rho=rec.rho_km_vals,
                                data_pow=rec.power_vals, verbose=False,
                                data_phase=rec.phase_rad_vals)

        grec = dc.DiffractionCorrection(gdata, 1.0, rng=rng)

        # Configurations for the plots.
        plt.rc('font', family='serif')
        plt.rc('font', size=10)
        plt.figure(figsize=(8.5, 11))

        # Title for the page.
        plt.suptitle("ModelFromGEO - Real Data - 200m Model", size=14)

        # Use gridspec to create a 4-by-1 group of plots.
        gs = gridspec.GridSpec(2, 1, hspace=0.0)

        plt.subplot(gs[0, 0])

        # Set up the tick parameters for the y axis.
        plt.tick_params(axis='y', which='both', left=True,
                        right=True, labelleft=True)
        plt.locator_params(axis='y', nbins=4)

        plt.title("Maxwell Ringlet - Rev007E Geometry Data")
        plt.tick_params(axis='x', which='both', bottom=False,
                        top=True, labelbottom=False)
        plt.locator_params(axis='x', nbins=8)

        # Set the label for the y-axis.
        plt.ylabel("Normalized Power")

        # Plot the modeled power.
        plt.plot(rec.rho_km_vals, rec.p_norm_vals, 'b', label="Actual Data")
        plt.plot(gdata.rho_km_vals, gdata.p_norm_vals, 'r',
                 label="Modeled Data")

        plt.xlim(rmin, rmax)
        plt.legend()

        plt.subplot(gs[1, 0])

        # Set up the tick parameters for the y axis.
        plt.tick_params(axis='y', which='both', left=True,
                        right=True, labelleft=True)
        plt.locator_params(axis='y', nbins=4)

        plt.xlabel("Ring Radius (km)")
        plt.ylabel("Normalized Power")
        plt.tick_params(axis='x', which='both', bottom=True,
                        top=False, labelbottom=True)
        plt.locator_params(axis='x', nbins=8)

        # Plot the reconstructed power.
        plt.plot(rec.rho_km_vals, rec.power_vals, 'r',
                 label="Real Reconstruction")
        plt.plot(grec.rho_km_vals, grec.power_vals, 'b',
                label="Reconstruction From Model")

        plt.xlim(rmin, rmax)
        plt.legend()

        # Append the plots to the pdf, close, and move on to the next one.
        pdf.savefig(bbox_inches="tight", pad_inches=1)
        plt.close()

if __name__ == "__main__":
    test_CompareTAU()
    test_ModelFromGEO()
