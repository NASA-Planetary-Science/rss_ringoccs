from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
sys.path.append("../")
from rss_ringoccs.tools.CSV_tools import ExtractCSVData
import rss_ringoccs.diffrec.advanced_tools as at
import rss_ringoccs.diffrec.diffraction_correction as dc


def test_ModelFromGEO():
    geo = "./Rev007E_X43_Maxwell_GEO.TAB"
    cal = "./Rev007E_X43_Maxwell_CAL.TAB"
    dlp = "./Rev007E_X43_Maxwell_DLP_500M.TAB"

    res = 0.2
    wav = 3.6e-5
    rho = 87510
    rng = [86000, 89000]
    dxs = 0.05
    outfile = "./figs/ModelFromGEO_Test.pdf"
    lbl = "Forward Model"
    lbf = "Fresnel Model"
    lba = "Actual Model"
    rlb = "Reconstruction from Forward Model"
    rfm = "Reconstruction from Fresnel Model"
    rmin = 87410
    rmax = 87610
    with PdfPages(outfile) as pdf:
        for model in ["Square Well", "Right Straight Edge",
                    "Left Straight Edge", "Delta Impulse"]:
            for opacity in {0.0, 0.3}:
                gdata = at.ModelFromGEO(geo, wav, res, rho, rng=rng,
                                        dx_km_desired=dxs, model=model,
                                        opacity=opacity)
                gdata_f = at.ModelFromGEO(geo, wav, res, rho, rng=rng,
                                          dx_km_desired=dxs, use_fresnel=True,
                                          model=model, opacity=opacity)
                rec = dc.DiffractionCorrection(gdata, res, rng=rng)
                rec_f = dc.DiffractionCorrection(gdata_f, res, rng=rng)

                # Configurations for the plots.
                plt.rc('font', family='serif')
                plt.rc('font', size=10)
                plt.figure(figsize=(8.5, 11))

                # Title for the page.
                plt.suptitle("ModelFromGEO - %s - 200m Model - Opacity: %.1f"
                             % (model, opacity), size=14)

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
                plt.plot(gdata_f.rho_km_vals, gdata_f.p_norm_vals, 'r', label=lbf)
                plt.plot(gdata.rho_km_vals, gdata.p_norm_vals, 'b', label=lbl)

                plt.xlim(rmin, rmax)
                plt.legend()
                if model == "Delta Impulse":
                    plt.ylim(0, 1.2)

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
                plt.plot(gdata.rho_km_vals, gdata.p_norm_actual_vals, 'g', label=lba)
                plt.plot(rec_f.rho_km_vals, rec_f.power_vals, 'r', label=rfm)
                plt.plot(rec.rho_km_vals, rec.power_vals, 'b', label=rlb)

                plt.xlim(rmin, rmax)
                plt.legend()

                # Append the plots to the pdf, close, and move on to the next one.
                pdf.savefig(bbox_inches="tight", pad_inches=1)
                plt.close()

        data = ExtractCSVData(geo, cal, dlp)
        rec = dc.DiffractionCorrection(data, 1.0, rng=[85000, 89000])
        gdata = at.ModelFromGEO(geo, wav, res, rho, rng=rng, dx_km_desired=dxs,
                                model="fromdata", data_rho=rec.rho_km_vals,
                                data_pow=rec.power_vals,
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
    test_ModelFromGEO()