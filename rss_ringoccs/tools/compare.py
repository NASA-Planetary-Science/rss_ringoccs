from rss_ringoccs import diffrec
from .CSV_tools import ExtractCSVData
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy import interpolate

def cringplots(rev, geo, cal, dlp, res, outfile="outfile.pdf",
               wtype="kbmd20", psitype="Fresnel4"):
    data = ExtractCSVData(geo, cal, dlp, verbose=False)
    TauInst = diffrec.DiffractionCorrection(data, res, wtype=wtype,
                                            rng=[74630, 90300], psitype=psitype)
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
            plt.suptitle("%s C-Ring Reconstructions: %s Resolution"
                         % (rev, res), size=14)
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
            plt.suptitle("%s C-Ring Reconstructions: %s Resolution"
                         % (rev, res), size=14)
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

def compare(NormDiff, geo, cal, dlp, tau, outfile, res=0.75, rng="all",
            wtype="kbmd20", norm=True, bfac=True, sigma=2.e-13, verbose=True,
            psitype="Fresnel8"):
    data = ExtractCSVData(geo, cal, dlp, tau=tau, verbose=verbose)
    rec = diffrec.DiffractionCorrection(NormDiff, res, rng=rng, wtype=wtype,
                                        fwd=True, norm=norm, bfac=bfac,
                                        sigma=sigma, verbose=verbose,
                                        psitype=psitype, write_file=False)
    rmin = np.min(rec.rho_km_vals)
    rmax = np.max(rec.rho_km_vals)
    nmin = np.min((rec.rho_km_vals >= rmin).nonzero())
    nmax = np.min((rec.rho_km_vals <= rmax).nonzero())

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
        p1 = rec.power_vals
        p2 = data.power_vals
        ymin = np.min([np.min(p1), np.min(p2)])*0.98
        ymax = np.max([np.max(p1), np.max(p2)])*1.02
        plt.xlim(rmin, rmax)
        plt.ylim(ymin, ymax)

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
    # Make sure rev, geo, cal, and dlp are strings.
    if (not isinstance(rev, str)):
        raise TypeError(
            "rev must be a string: 'Rev###'\n"
            "\tYour input has type: %s\n"
            "\tInput should have type: str\n"
            % (type(rev).__name__)
        )
    elif (not isinstance(geo, str)):
        raise TypeError(
            "geo must be a string: '/path/to/geo'\n"
            "\tYour input has type: %s\n"
            "\tInput should have type: str\n"
            % (type(geo).__name__)
        )
    elif (not isinstance(cal, str)):
        raise TypeError(
            "cal must be a string: '/path/to/cal'\n"
            "\tYour input has type: %s\n"
            "\tInput should have type: str\n"
            % (type(cal).__name__)
        )
    elif (not isinstance(dlp, str)):
        raise TypeError(
            "dlp must be a string: '/path/to/dlp'\n"
            "\tYour input has type: %s\n"
            "\tInput should have type: str\n"
            % (type(dlp).__name__)
        )
    elif (not isinstance(outfile, str)):
        raise TypeError(
            "outfile must be a string: '/path/to/outfile'\n"
            "\tYour input has type: %s\n"
            "\tInput should have type: str\n"
            % (type(outfile).__name__)
        )
    elif (not isinstance(ymin, float)):
        try:
            ymin = float(ymin)
        except:
            raise TypeError(
                "ymin must be a floating point number.\n"
                "\tYour input has type: %s\n"
                "\tInput should have type: float\n"
                % (type(ymin).__name__)
            )
    else:
        pass

    if tau:
        if (not isinstance(tau, str)):
            raise TypeError(
                "tau must be a string: '/path/to/tau'\n"
                "\tYour input has type: %s\n"
                "\tInput should have type: str\n"
                % (type(tau).__name__)
            )
        else:
            pass
    else:
        pass

    # Make sure that verbose is a boolean.
    if not isinstance(verbose, bool):
        raise TypeError(
            "\n\tverbose must be Boolean: True/False\n"
            "\tYour input has type: %s\n"
            "\tInput should have type: bool\n"
            "\tSet verbose=True or verbose=False\n"
            % (type(verbose).__name__)
        )
    else:
        pass

    # Check that the input resolution is a positive floating point number.
    if (not isinstance(res, list)):
        try:
            res = float(res)
            res = [res]
        except (TypeError, ValueError):
            raise TypeError(
                "\n\tres must be a list of positive floating point numbers\n"
                "\tYour input has type: %s\n"
                "\tInput should have type: float\n"
                % (type(res).__name__)
            )
    else:
        pass

    for x in res:
        if (x <= 0.0):
            raise ValueError(
                "\n\tres must be a list of positive floating point numbers\n"
                "\tYour requested resolutions (km): %f\n"
                % (res)
            )
        else:
            pass

    # Ensure that the normalize boolean has a legal value.
    if not isinstance(norm, bool):
        raise TypeError(
            "\n\tnorm must be Boolean: True/False\n"
            "\tYour input has type: %s\n"
            "\tInput should have type: bool\n"
            "\tSet norm=True or norm=False\n"
            % (type(norm).__name__)
        )
    else:
        pass

    # Make sure that bfac is a boolean.
    if not isinstance(bfac, bool):
        raise TypeError(
            "\n\tbfac must be Boolean: True/False\n"
            "\tYour input has type: %s\n"
            "\tInput should have type: bool\n"
            "\tSet bfac=True or bfac=False\n"
            % (type(bfac).__name__)
        )
    else:
        pass
        
    # Check that res_factor is a floating point number.
    if (not isinstance(res_factor, float)):
        try:
            res_factor = float(res_factor)
        except (TypeError, ValueError):
            raise TypeError(
                "\n\tres_factor must be a positive floating point number.\n"
                "\tYour input has type: %s\n"
                % (type(res_factor).__name__)
            )
    else:
        pass

    if (res_factor <= 0.0):
            raise ValueError(
                "\n\tres_factor must be a positive floating point number.\n"
                "\tYour input: %f\n"
                % (res_factor)
            )
    else:
        pass

    # Check that the Allen Deviation is a legal value.
    if (not isinstance(sigma, float)):
        try:
            sigma = float(sigma)
        except (TypeError, ValueError):
            raise TypeError(
                "\n\tsigma must be a positive floating point number.\n"
                "\tYour input has type: %s\n"
                % (type(sigma).__name__)
            )
    else:
            pass

    if (np.min(sigma) <= 0.0):
            raise ValueError(
                "\n\tsigma must be a positive floating point number.\n"
                "\tYour input: %f\n"
                % (sigma)
            )
    else:
        pass


    data = ExtractCSVData(geo, cal, dlp, tau=tau, verbose=False)
    N_Plots = len(res)
    N_Pages = int(N_Plots/4.0)

    with PdfPages(outfile) as pdf:
        for i_page in range(N_Pages):
            plt.rc('font', family='serif')
            plt.rc('font', size=10)
            plt.figure(figsize=(8.5, 11))
            plt.suptitle("Resolution Comparison: %s" % (rev), size=14)
            gs = gridspec.GridSpec(4, 1, hspace=0.0)

            # Perform Reconstructions
            for i in range(4):
                i_res = int(4*i_page+i)
                sres = str(res[i_res])+"km Reconstruction"
                rec = diffrec.DiffractionCorrection(data, res[i_res],
                                                    wtype=wtype, rng=rng,
                                                    psitype=psitype, norm=norm,
                                                    bfac=bfac, sigma=sigma,
                                                    res_factor=res_factor,
                                                    verbose=verbose)
                plt.subplot(gs[i, 0])
                plt.tick_params(axis='y', which='both', left=True,
                                right=True, labelleft=True)
                plt.locator_params(axis='y', nbins=4)
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
                plt.ylabel("Normalized Power")
                plt.plot(rec.rho_km_vals, rec.power_vals, 'b', label=sres)
                plt.plot(data.tau_rho, data.power_vals, 'r', label="PDS")
                rmin = np.min(rec.rho_km_vals)
                rmax = np.max(rec.rho_km_vals)
                plt.xlim(rmin, rmax)
                plt.ylim(ymin, ymax)
                plt.legend()

            pdf.savefig(bbox_inches="tight", pad_inches=1)
            plt.close()

        if ((N_Plots % 4) != 0):
            plt.rc('font', family='serif')
            plt.rc('font', size=10)
            plt.figure(figsize=(8.5, 11))
            plt.suptitle("Resolution Comparison: %s" % (rev), size=14)
            gs = gridspec.GridSpec(4, 1, hspace=0.0)

            for i in range(int(N_Plots % 4)):
                i_res = int(4*N_Pages+i)
                sres = str(res[i_res])+"km Reconstruction"
                rec = diffrec.DiffractionCorrection(data, res[i_res],
                                                    wtype=wtype, rng=rng,
                                                    psitype=psitype, norm=norm,
                                                    bfac=bfac, sigma=sigma,
                                                    res_factor=res_factor,
                                                    verbose=verbose)
                plt.subplot(gs[i, 0])
                plt.tick_params(axis='y', which='both', left=True,
                                right=True, labelleft=True)
                plt.locator_params(axis='y', nbins=4)
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
                plt.ylabel("Normalized Power")
                plt.plot(rec.rho_km_vals, rec.power_vals, 'b', label=sres)
                plt.plot(data.tau_rho, data.power_vals, 'r', label="PDS")
                rmin = np.min(rec.rho_km_vals)
                rmax = np.max(rec.rho_km_vals)
                plt.xlim(rmin, rmax)
                plt.ylim(ymin, ymax)
                plt.legend()
            pdf.savefig(bbox_inches="tight", pad_inches=1)
            plt.close()
        else:
            pass
