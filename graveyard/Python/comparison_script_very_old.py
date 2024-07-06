"""
################################################################################
#                                   LICENSE                                    #
################################################################################
#   This file is part of rss_ringoccs.                                         #
#                                                                              #
#   rss_ringoccs is free software: you can redistribute it and/or              #
#   modify it under the terms of the GNU General Public License as published   #
#   by the Free Software Foundation, either version 3 of the License, or       #
#   (at your option) any later version.                                        #
#                                                                              #
#   rss_ringoccs is distributed in the hope that it will be useful             #
#   but WITHOUT ANY WARRANTY# without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#   GNU General Public License for more details.                               #
#                                                                              #
#   You should have received a copy of the GNU General Public License          #
#   along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.     #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
# pylint: disable = invalid-name
# pylint: disable = consider-using-f-string
# pylint: disable = too-many-locals
# pylint: disable = too-many-statements
import numpy
from scipy import interpolate
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib import gridspec
import spiceypy
from .date_string import date_string

def make_plots(c1, c2, title):
    """
        Makes comparison plots for Cassini data.
    """
    date = date_string()
    savestring = date + title.replace(" ", "_") + ".pdf"

    rho1 = c1.rho_km_vals
    rho2 = c2.rho_km_vals

    dx1 = rho1[1]-rho1[0]
    dx2 = rho2[1]-rho2[0]

    rmin = numpy.max([numpy.min(rho1), numpy.min(rho2)])
    rmax = numpy.min([numpy.max(rho1), numpy.max(rho2)])

    dx = numpy.min([dx1, dx2])
    n = int(numpy.floor((rmax - rmin) / dx + 1.0))
    r = rmin + numpy.arange(n) * dx

    power1 = c1.power_vals
    power2 = c2.power_vals
    phase1 = c1.phase_vals
    phase2 = c2.phase_vals
    tau1 = c1.tau_vals
    tau2 = c2.tau_vals

    epower1 = c1.tau_power
    epower2 = c2.tau_power
    ephase1 = c1.tau_phase
    ephase2 = c2.tau_phase
    etau1 = c1.tau_tau
    etau2 = c2.tau_tau

    p_interp1 = interpolate.interp1d(rho1,power1,kind='linear')
    phase_interp1 = interpolate.interp1d(rho1,phase1,kind='linear')
    tau_interp1 = interpolate.interp1d(rho1,tau1,kind='linear')
    p_interp2 = interpolate.interp1d(rho2,power2,kind='linear')
    phase_interp2 = interpolate.interp1d(rho2,phase2,kind='linear')
    tau_interp2 = interpolate.interp1d(rho2,tau2,kind='linear')

    wtype1 = c1.wtype
    wtype2 = c2.wtype

    res1 = c1.res
    res2 = c2.res

    power1 = p_interp1(r)
    power2 = p_interp2(r)

    phase1 = phase_interp1(r)
    phase2 = phase_interp2(r)

    tau1 = tau_interp1(r)
    tau2 = tau_interp2(r)

    del p_interp1, p_interp2, phase_interp1
    del phase_interp2, tau_interp1, tau_interp2

    p_interp1 = interpolate.interp1d(rho1,epower1,kind='linear')
    phase_interp1 = interpolate.interp1d(rho1,ephase1,kind='linear')
    tau_interp1 = interpolate.interp1d(rho1,etau1,kind='linear')
    p_interp2 = interpolate.interp1d(rho2,epower2,kind='linear')
    phase_interp2 = interpolate.interp1d(rho2,ephase2,kind='linear')
    tau_interp2 = interpolate.interp1d(rho2,etau2,kind='linear')

    epower1 = p_interp1(r)
    epower2 = p_interp2(r)
    ephase1 = phase_interp1(r)
    ephase2 = phase_interp2(r)
    etau1 = tau_interp1(r)
    etau2 = tau_interp2(r)

    phase1 = phase1*spiceypy.dpr()
    phase2 = phase2*spiceypy.dpr()
    ephase1 = ephase1*spiceypy.dpr()
    ephase2 = ephase2*spiceypy.dpr()

    del p_interp1, p_interp2, phase_interp1
    del phase_interp2, tau_interp1, tau_interp2

    with PdfPages(savestring) as pdf:
        plt.rc('text', usetex = True)
        plt.rc('font', family = 'serif')
        plt.rc('font', size = 10)
        plt.figure(figsize = (8.5, 11))
        plt.suptitle(title.replace("_", r"\_"), size = 14)
        gs = gridspec.GridSpec(3, 1, hspace = 0.33)
        p1 = gridspec.GridSpecFromSubplotSpec(
            2, 2, gs[0, 0], hspace = 0.0, wspace = 0.0
        )

        p2 = gridspec.GridSpecFromSubplotSpec(
            2, 2, gs[1, 0], hspace = 0.0, wspace = 0.0
        )

        p3 = gridspec.GridSpecFromSubplotSpec(
            2, 2, gs[2, 0], hspace = 0.0, wspace = 0.0
        )

        xmin, xmax = numpy.min(r), numpy.max(r)
        xlim   = (xmin, xmax)
        ymin1a = numpy.min([power1, power2, epower1, epower2]) * 0.80
        ymax1a = numpy.max([power1, power2, epower1, epower2]) * 1.20

        ymin1b = numpy.max(
            [numpy.min(power1 - epower1), numpy.min(power2 - epower2)]
        )

        ymax1b = numpy.min(
            [numpy.max(power1 - epower1), numpy.max(power2 - epower2)]
        )

        ythin1 = numpy.max(numpy.abs([ymin1b,ymax1b])) * 2.0
        ylim1a = (ymin1a, ymax1a)
        ylim1b = (-ythin1, ythin1)

        del ymin1a, ymax1a, ymin1b, ymax1b

        ymin2a = numpy.min([phase1, phase2, ephase1, ephase2])
        ymin2b = numpy.min([phase1 - ephase1, phase2 - ephase2])
        ymax2a = numpy.max([phase1, phase2, ephase1, ephase2])
        ymax2b = numpy.max([phase1 - ephase1, phase2 - ephase2])

        ylim2a = (ymin2a, ymax2a)
        ylim2b = (ymin2b, ymax2b)

        del ymin2a, ymax2a, ymin2b, ymax2b

        ymin3b = numpy.max([numpy.min(tau1 - etau1), numpy.min(tau2 - etau2)])
        ymax3a = numpy.max([tau1, tau2, etau1, etau2])
        ymax3b = numpy.min([numpy.max(tau1 - etau1), numpy.max(tau2 - etau2)])
        ythin2 = numpy.max(numpy.abs([ymin3b, ymax3b])) * 2.0
        ylim3a = (0.0, ymax3a)
        ylim3b = (-ythin2, ythin2)

        del ymin3b, ymax3a, ymax3b

        plt.subplot(p1[0.0, 0.0])
        plt.title("Power - %s - %.3f km" % (wtype1, res1))
        plt.ylabel("Power")
        plt.xlim(xlim)
        plt.ylim(ylim1a)
        plt.tick_params(
            axis = 'y', which = 'both',
            left = True, right = False, labelleft = True
        )

        plt.tick_params(
            axis = 'x', which = 'both',
            bottom = False, top = True, labelbottom = False
        )

        plt.yticks(numpy.arange(0.0, max(ylim1a), 0.25))

        plt.plot(r, power1, label = "TC", linewidth = 1)
        plt.plot(r, epower1, label = "PDS", linewidth = 1)

        plt.legend()
        plt.subplot(p1[1, 0])
        plt.tick_params(
            axis = 'y', which = 'both',
            left = True, right = False, labelleft = True
        )

        plt.tick_params(
            axis = 'x', which = 'both',
            bottom = True, top = False, labelbottom = True
        )

        plt.yticks(numpy.arange(-ythin1, 0.51*ythin1, 0.5*ythin1))
        plt.xlabel('Ring Radius (km)')
        plt.ylabel("Difference")
        plt.xlim(xlim)
        plt.ylim(ylim1b)
        plt.plot(r, power1 - epower1, linewidth = 1)
        plt.subplot(p1[0, 1])
        plt.title("Power - %s - %.3f km" % (wtype2, res2))

        plt.tick_params(
            axis = 'y', which = 'both',
            left = False, right = True, labelleft = False
        )

        plt.tick_params(
            axis = 'x', which = 'both',
            bottom = False, top = True, labelbottom = False
        )

        plt.yticks(numpy.arange(0, max(ylim1a), 0.25))
        plt.xlim(xlim)
        plt.ylim(ylim1a)

        plt.plot(r, power2, label = "TC", linewidth = 1)
        plt.plot(r, epower2, label = "PDS", linewidth = 1)

        plt.legend()
        plt.subplot(p1[1, 1])
        plt.yticks(numpy.arange(-ythin1, 0.51*ythin1, 0.5*ythin1))

        plt.tick_params(
            axis = 'y', which = 'both',
            left = False, right = True, labelleft = False
        )

        plt.tick_params(
            axis = 'x', which = 'both',
            bottom = True, top = False, labelbottom = True
        )

        plt.xlabel('Ring Radius (km)')
        plt.xlim(xlim)
        plt.ylim(ylim1b)
        plt.plot(r, power2 - epower2, linewidth = 1)

        plt.subplot(p2[0,0])
        plt.title("Phase - %s - %.3f km" % (wtype1,res1))
        plt.ylabel("Phase (Degrees)")
        plt.tick_params(
            axis = 'x', which = 'both',
            bottom = False, top = True, labelbottom = False
        )

        plt.xlim(xlim)
        plt.ylim(ylim2a)

        plt.plot(r, phase1, label = "TC", linewidth = 1)
        plt.plot(r, ephase1, label = "PDS", linewidth = 1)

        plt.legend()
        plt.subplot(p2[1, 0])
        plt.xlabel('Ring Radius (km)')
        plt.ylabel("Difference")
        plt.xlim(xlim)
        plt.ylim(ylim2b)

        plt.plot(r, phase1 - ephase1, linewidth = 1)
        plt.subplot(p2[0,1])

        plt.title("Phase - %s - %.3f km" % (wtype2, res2))

        plt.tick_params(
            axis = 'y', which = 'both',
            left = False, right = True, labelleft = False
        )

        plt.tick_params(
            axis = 'x', which = 'both',
            bottom = False, top = True, labelbottom = False
        )

        plt.xlim(xlim)
        plt.ylim(ylim2a)
        plt.plot(r, phase2, label = "TC", linewidth = 1)
        plt.plot(r, ephase2, label = "PDS", linewidth = 1)
        plt.legend()
        plt.subplot(p2[1,1])
        plt.xlabel('Ring Radius (km)')

        plt.tick_params(
            axis = 'y', which = 'both',
            left = False, right = True, labelleft = False
        )

        plt.xlim(xlim)
        plt.ylim(ylim2b)
        plt.plot(r, phase2 - ephase2, linewidth = 1)

        plt.subplot(p3[0,0])
        plt.title("Optical Depth - %s - %.3f km" % (wtype1,res1))
        plt.ylabel("Optical Depth")
        plt.tick_params(
            axis = 'x', which = 'both',
            bottom = False, top = True, labelbottom = False
        )

        plt.tick_params(
            axis = 'y', which = 'both',
            left = True, right = False, labelleft = True
        )

        plt.yticks(numpy.arange(0, max(ylim3a), 0.5))
        plt.xlim(xlim)
        plt.ylim(ylim3a)
        plt.plot(r, tau1, label = "TC", linewidth = 1)
        plt.plot(r, etau1, label = "PDS", linewidth = 1)
        plt.legend()
        plt.subplot(p3[1, 0])
        plt.xlabel("Ring Radius (km)")
        plt.ylabel("Optical Depth")
        plt.xlim(xlim)
        plt.ylim(ylim3b)

        plt.tick_params(
            axis = 'y', which = 'both',
            left = True, right = False, labelleft = True
        )

        plt.yticks(numpy.arange(-ythin2, 0.51*ythin2, 0.5*ythin2))
        plt.plot(r, tau1 - etau1, linewidth = 1)
        plt.subplot(p3[0, 1])

        plt.title("Optical Depth - %s - %.3f km" % (wtype2, res2))

        plt.tick_params(
            axis = 'y', which = 'both',
            left = False, right = True, labelleft = False
        )

        plt.tick_params(
            axis = 'x', which = 'both',
            bottom = False, top = True, labelbottom = False
        )

        plt.yticks(numpy.arange(0, max(ylim3a), 0.5))
        plt.xlim(xlim)
        plt.ylim(ylim3a)
        plt.plot(r, tau2, label = "TC", linewidth = 1)
        plt.plot(r, etau2, label = "PDS", linewidth = 1)
        plt.legend()
        plt.subplot(p3[1,1])
        plt.xlabel('Ring Radius (km)')

        plt.tick_params(
            axis = 'y', which = 'both',
            left = False, right = True, labelleft = False
        )

        plt.yticks(numpy.arange(-ythin2, 0.51*ythin2, 0.5*ythin2))
        plt.xlim(xlim)
        plt.ylim(ylim3b)
        plt.plot(r, tau2 - etau2, linewidth = 1)

        # Saves the current figure into a pdf page.
        pdf.savefig(bbox_inches = "tight", pad_inches = 1)
        plt.close()
        