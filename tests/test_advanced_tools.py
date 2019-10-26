from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
sys.path.append("../")
import rss_ringoccs.diffrec.advanced_tools as at
import rss_ringoccs.diffrec.diffraction_correction as dc


if __name__ == "__main__":
    geo = "./Rev007E_X43_Maxwell_GEO.TAB"

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
            gdata = at.ModelFromGEO(geo, wav, res, rho, rng=rng,
                                    dx_km_desired=dxs, model=model,
                                    verbose=False)
            gdata_f = at.ModelFromGEO(geo, wav, res, rho, rng=rng,
                                    dx_km_desired=dxs, use_fresnel=True,
                                    model=model,verbose=False)
            rec = dc.DiffractionCorrection(gdata,res,rng=rng,verbose=False)
            rec_f = dc.DiffractionCorrection(gdata_f,res,rng=rng,verbose=False)

            # Configurations for the plots.
            plt.rc('font', family='serif')
            plt.rc('font', size=10)
            plt.figure(figsize=(8.5, 11))

            # Title for the page.
            plt.suptitle("%s - 200m Model" % model, size=14)

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

            # Plot the reconstructed power.
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