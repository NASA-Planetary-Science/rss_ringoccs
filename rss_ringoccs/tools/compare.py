def compare(NormDiff, geo, cal, dlp, tau, outfile, res=0.75, rng="all",
            wtype="kbmd20", norm=True, bfac=True, sigma=2.e-13, verbose=True,
            psitype="Fresnel8"):

    # Retrieve the data from the input CSV files.
    data = CSV_tools.ExtractCSVData(geo, cal, dlp, tau=tau, verbose=verbose)

    # Perform diffraction reconstruction on the raw data.
    rec = diffrec.DiffractionCorrection(NormDiff, res, rng=rng, wtype=wtype,
                                        fwd=True, norm=norm, bfac=bfac,
                                        sigma=sigma, verbose=verbose,
                                        psitype=psitype, write_file=False)

    # Compute the range for the x-axis.
    rmin = np.min(rec.rho_km_vals)
    rmax = np.max(rec.rho_km_vals)

    # Compute the indices corresponding to the range of the x-axis.
    nmin = np.min((rec.rho_km_vals >= rmin).nonzero())
    nmax = np.min((rec.rho_km_vals <= rmax).nonzero())

    with PdfPages(outfile) as pdf:
        # Configuration for the plots.
        plt.rc('font', family='serif')
        plt.rc('font', size=10)
        plt.figure(figsize=(8.5, 11))

        # Title for the page.
        plt.suptitle("TAU Parameters", size=14)

        # Use gridspec to create a 3-by-2 group of plots.
        gs = gridspec.GridSpec(3, 2, wspace=0.0, hspace=0.0)

        # Subplot for Normalized Power.
        plt.subplot(gs[0, 0])

        # Set tick parameters for the x and y axes.
        plt.tick_params(axis='y', which='both', left=True,
                        right=False, labelleft=True)
        plt.tick_params(axis='x', which='both', bottom=False,
                        top=True, labelbottom=False)

        # Add the title for the plots and labels.
        plt.title("Comparison Plots")
        plt.ylabel('Normalized Power')

        # Plot reconstructed power and the power being compared against.
        plt.plot(rec.rho_km_vals, rec.power_vals, 'b')
        plt.plot(data.tau_rho, data.power_vals, 'g')
        p1 = rec.power_vals
        p2 = data.power_vals

        # Compute the range of the y-axis.
        ymin = np.min([np.min(p1), np.min(p2)])*0.98
        ymax = np.max([np.max(p1), np.max(p2)])*1.02
        plt.xlim(rmin, rmax)
        plt.ylim(ymin, ymax)

        # Plot Difference in Normalized Power.
        plt.subplot(gs[0, 1])
        interp = interpolate.interp1d(data.tau_rho, data.power_vals)
        diff = rec.power_vals - interp(rec.rho_km_vals)

        # Set tick parameters for the difference plot.
        plt.tick_params(axis='y', which='both', left=False,
                        right=True, labelleft=False, labelright=True)
        plt.tick_params(axis='x', which='both', bottom=False,
                        top=True, labelbottom=False)

        # Add a title, and make the plots.
        plt.title("Difference Plots")
        plt.plot(rec.rho_km_vals, diff, 'r')

        # Set the range for the x and y axes.
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
        p1 = rec.tau_vals
        p2 = data.tau_vals
        ymin = np.min([np.min(p1), np.min(p2)])*0.98
        ymax = np.max([np.max(p1), np.max(p2)])*1.02
        plt.xlim(rmin, rmax)
        plt.ylim(ymin, ymax)

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
        p1 = rec.phase_vals
        p2 = data.phase_vals
        ymin = np.min([np.min(p1), np.min(p2)])*0.98
        ymax = np.max([np.max(p1), np.max(p2)])*1.02
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
        plt.ylim(np.min(diff), np.max(diff))

        pdf.savefig(bbox_inches="tight", pad_inches=1)
        plt.close()

def galleryplots(rev, geo, cal, dlp, tau=None, res=[1.0], rng="all",
                 wtype="kbmd20", psitype="Fresnel4", norm=True, bfac=True,
                 sigma=2.e-13, verbose=True, res_factor=0.75,
                 outfile="galleryplot.pdf", ymin=-0.2, ymax=1.4):
    """
        Purpose:
            Create a set of plots of the same ring feature at
            various resolutions as specified by the user.
        Arguments:
            :rev (*str*):
                The rev number. Ex: rev = "Rev007"
            :geo (*str*):
                The location of the geo file.
                Ex: geo = "/path/to/geo"
            :cal (*str*)
                The location of the cal file.
                Ex: dlp = "/path/to/cal"
            :dlp (*str*):
                The location of the dlp file.
                Ex: dlp = "/path/to/dlp"
        Keywords:
            :tau (*str*):
                The location of the tau file. If set,
                this plots the PDS power, as well as the
                user reconstructed power.
                Ex: tau = "/path/to/tau"
            :res (*list*):
                The set of requested resolution to process
                and plot. The values should be floating
                point values and take the sampling theorem
                into consideration.
                Ex: res = [0.5, 0.7, 1.0, 1.2]
            :rng (*list* or *str*):
                The requested range for diffraction correction.
                Preferred input is rng = [a,b]. Arrays are
                allowed and the range will be set as:

                |    rng = [MIN(array), MAX(array)]

                Finally, certain strings containing a few of the
                regions of interests within the rings of Saturn
                are allowed. Permissable strings are:

                |    'all'             [1.0, 400000.0]
                |    'cringripples'    [77690.0, 77760.0]
                |    'encke'           [132900.0, 134200.0]
                |    'enckegap'        [132900.0, 134200.0]
                |    'janusepimetheus' [96200.0, 96800.0]
                |    'maxwell'         [87410.0, 87610.0]
                |    'maxwellringlet'  [87410.0, 87610.0]
                |    'titan'           [77870.0, 77930.0]
                |    'titanringlet'    [77870.0, 77930.0]
                |    'huygens'         [117650.0, 117950.0]
                |    'huygensringlet'  [117650.0, 117950.0]

                Strings are neither case nor space sensitive.
                For other planets use rng = [a,b]. Default value
                is set to 'all' which processes [1, 400000]
                Values MUST be set in kilometers.
            :wtype (*str):
                The requested tapering function for diffraction
                correction. A string with several allowed inputs:

                |    'rect'      Rectangular Window.
                |    'coss'      Squared Cosine Window.
                |    'kb20'      Kaiser-Bessel 2.0 Window.
                |    'kb25'      Kaiser-Bessel 2.5 Window.
                |    'kb35'      Kaiser-Bessel 3.5 Window.
                |    'kbmd20'    Modified kb20 Window.
                |    'kbmd25'    Modified kb25 Window.

                The variable is neither case nor space sensitive.
                Default window is set to 'kb25'. See window_functions
                submodule for further documentation.
            :norm (*bool*):
                A Boolean for determining whether or not the
                reconstructed complex transmittance is normalize
                by the window width. This normalization is the
                complex transmittance that is computed by using
                free space divided by the complex transmittance
                that is computed using free space weighted by the
                selected tapering function. Default is True.
            :bfac (*bool*):
                A Boolean for determining whether or not the
                'b' factor in the window width computation is
                used. This is equivalent to setting the Allen
                Deviation for the spacecraft to a positive value
                or to zero. If set to False, the Allen Deviation
                is assumed to be zero. If set to True the Allen
                Deviation is set to 2e-13, or whichever number you
                wish to specify in the sigma keyword (See below).
                Default is True.
            :sigma (*float*):
                The Allen deviation for the spacecraft. If the bfac
                keyword (See above) is set to False, this is ignored.
                If bfac is set to True, and sigma is NOT specified,
                then sigma=2e-13 will be used, which is the Allen
                deviation for Cassini with 1 second integration time.
                For spacecraft other than Cassini, you should provide
                the Allen deviation yourself. Default is sigma=2e-13
            :psitype (*str*):
                A string for determining what approximation to the
                geometrical 'psi' function is used. Several strings
                are allowed:

                |    'full'      No Approximation is applied.
                |    'MTR2'      Second Order Series from MTR86.
                |    'MTR3'      Third Order Series from MTR86.
                |    'MTR4'      Fourth Order Series from MTR86.
                |    'Fresnel'   Standard Fresnel approximation.

                The variable is neither case nor space sensitive.
                Default is set to 'full'.
            :verbose (*bool*):
                A Boolean for determining if various pieces of
                information are printed to the screen or not.
                Default is False.
            :outfile (*str*):
                Path to the output folder and the name of the
                pdf that is to be created.
                Ex: outfile = "/path/to/outfile.pdf"
            :res_factor (*float*):
                Floating point number used as a scale factor
                for the resolution for the sake of consistency
                with the PDS results. The definition of
                resolution adopted in the PDS and the
                definition specified in MTR86 differs by
                a scale of about 0.75. To skip this, set
                res_factor = 1.0.
            :ymin (*float*):
                The minimum y value to be plotted.
                Ex: ymin = -0.2
            :ymax (*float*):
                The maximum y value to be plotted.
                Ex: ymax = 1.5
    """
    # Check that the input variables contain legal values.
    fname = "tools.compare.galleryplots"
    error_check.check_type(rev, str, "rev", fname)
    error_check.check_type(geo, str, "geo", fname)
    error_check.check_type(cal, str, "cal", fname)
    error_check.check_type(dlp, str, "dlp", fname)
    error_check.check_type(norm, bool, "norm", fname)
    error_check.check_type(bfac, bool, "bfac", fname)
    error_check.check_type(outfile, str, "outfile", fname)
    error_check.check_type(verbose, bool, "verbose", fname)
    ymin = error_check.check_type_and_convert(ymin, float, "ymin", fname)
    ymax = error_check.check_type_and_convert(ymax, float, "ymax", fname)
    sigma = error_check.check_type_and_convert(sigma, float, "sigma", fname)
    res_factor = error_check.check_type_and_convert(res_factor, float,
                                                    "res_factor", fname)
    
    error_check.check_positive(sigma, "sigma", fname)
    error_check.check_positive(res_factor, "res_factor", fname)

    psitype = error_check.check_psitype(psitype, fname)

    # Check that the input resolution is a positive floating point number.
    if (not isinstance(res, list)):
        res = error_check.check_type_and_convert(res, float, "res", fname)
        res = [res]
    else:
        pass

    for x in res:
        error_check.check_positive(x, "res", fname)
    
    if tau:
        error_check.check_type(tau, str, "tau", fname)

    # Extract the data from the given CSV files.
    data = CSV_tools.ExtractCSVData(geo, cal, dlp, tau=tau, verbose=verbose)

    # Total number of plots, and total number of pages.
    N_Plots = len(res)
    N_Pages = int(N_Plots/4.0)

    with PdfPages(outfile) as pdf:
        for i_page in range(N_Pages):
            # Configurations for the plots.
            plt.rc('font', family='serif')
            plt.rc('font', size=10)
            plt.figure(figsize=(8.5, 11))

            # Title for the page.
            plt.suptitle("Resolution Comparison: %s" % (rev), size=14)

            # Use gridspec to create a 4-by-1 group of plots.
            gs = gridspec.GridSpec(4, 1, hspace=0.0)

            # Perform Reconstructions
            for i in range(4):
                i_res = int(4*i_page+i)
                sres = str(res[i_res])+"km Reconstruction"

                # Use the DiffractionCorrection class to get the output data.
                rec = diffrec.DiffractionCorrection(data, res[i_res],
                                                    wtype=wtype, rng=rng,
                                                    psitype=psitype, norm=norm,
                                                    bfac=bfac, sigma=sigma,
                                                    res_factor=res_factor,
                                                    verbose=verbose)
                plt.subplot(gs[i, 0])

                # Set up the tick parameters for the y axis.
                plt.tick_params(axis='y', which='both', left=True,
                                right=True, labelleft=True)
                plt.locator_params(axis='y', nbins=4)

                # x-axis ticks depend on if this is the top, middle, or bottom.
                if (i == 0):
                    start = str(np.min(rec.rho_km_vals))
                    end = str(np.max(rec.rho_km_vals))
                    plt.title("Comparison Plots: %skm to %skm" % (start, end))
                    plt.tick_params(axis='x', which='both', bottom=False,
                                    top=True, labelbottom=False)
                    plt.locator_params(axis='x', nbins=8)
                elif (i == 3):
                    plt.xlabel("Ring Radius (km)")
                    plt.tick_params(axis='x', which='both', bottom=True,
                                    top=False, labelbottom=True)
                    plt.locator_params(axis='x', nbins=8)
                else:
                    plt.tick_params(axis='x', which='both', bottom=False,
                                    top=False, labelbottom=False)

                # Set the label for the y-axis.
                plt.ylabel("Normalized Power")

                # Plot the reconstructed power.
                plt.plot(rec.rho_km_vals, rec.power_vals, 'b', label=sres)

                # If a tau file is provided, plot the results it contains.
                if tau:
                    plt.plot(data.tau_rho, data.power_vals, 'r', label="PDS")

                # Set the x and y range for the plots.
                rmin = np.min(rec.rho_km_vals)
                rmax = np.max(rec.rho_km_vals)
                plt.xlim(rmin, rmax)
                plt.ylim(ymin, ymax)
                plt.legend()

            # Append the plots to the pdf, close, and move on to the next one.
            pdf.savefig(bbox_inches="tight", pad_inches=1)
            plt.close()

        if ((N_Plots % 4) != 0):
            # Configurations for the plots.
            plt.rc('font', family='serif')
            plt.rc('font', size=10)
            plt.figure(figsize=(8.5, 11))

            # Title for the page.
            plt.suptitle("Resolution Comparison: %s" % (rev), size=14)

            # Use gridspec to create a 4-by-1 group of plots.
            gs = gridspec.GridSpec(4, 1, hspace=0.0)

            for i in range(int(N_Plots % 4)):
                i_res = int(4*N_Pages+i)
                sres = str(res[i_res])+"km Reconstruction"

                # Perform diffraction correction on the last sets.
                rec = diffrec.DiffractionCorrection(data, res[i_res],
                                                    wtype=wtype, rng=rng,
                                                    psitype=psitype, norm=norm,
                                                    bfac=bfac, sigma=sigma,
                                                    res_factor=res_factor,
                                                    verbose=verbose)
                plt.subplot(gs[i, 0])

                # Set tick parameters for the y axis.
                plt.tick_params(axis='y', which='both', left=True,
                                right=True, labelleft=True)
                plt.locator_params(axis='y', nbins=4)

                # Set tick parameters for the x-axis.
                if (i == 0):
                    start = str(np.min(rec.rho_km_vals))
                    end = str(np.max(rec.rho_km_vals))
                    plt.title("Comparison Plots: %skm to %skm" % (start, end))
                    plt.tick_params(axis='x', which='both', bottom=False,
                                    top=True, labelbottom=False)
                    plt.locator_params(axis='x', nbins=8)
                elif (i == (N_Plots % 4) - 1):
                    plt.xlabel("Ring Radius (km)")
                    plt.tick_params(axis='x', which='both', bottom=True,
                                    top=False, labelbottom=True)
                    plt.locator_params(axis='x', nbins=8)
                else:
                    plt.tick_params(axis='x', which='both', bottom=False,
                                    top=False, labelbottom=False)
                
                # Add labels to the y-axis.
                plt.ylabel("Normalized Power")

                # Plot the results.
                plt.plot(rec.rho_km_vals, rec.power_vals, 'b', label=sres)

                # If a tau file is provided, plot the results it contains.
                if tau:
                    plt.plot(data.tau_rho, data.power_vals, 'r', label="PDS")
                
                # Set the range for the x and y axes.
                rmin = np.min(rec.rho_km_vals)
                rmax = np.max(rec.rho_km_vals)
                plt.xlim(rmin, rmax)
                plt.ylim(ymin, ymax)
                plt.legend()
            
            # Save the plots and close the PDF.
            pdf.savefig(bbox_inches="tight", pad_inches=1)
            plt.close()
        else:
            pass
