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
from rss_ringoccs import DiffractionCorrection, ExtractCSVData
import numpy
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from .date_string import date_string

def rev007_window():
    """
        Gets the data for the plots.
    """
    geo007 = 'Rev007E_RSS_2005_123_X43_E/RSS_2005_123_X43_E_GEO.TAB'
    cal007 = 'Rev007E_RSS_2005_123_X43_E/RSS_2005_123_X43_E_CAL.TAB'
    dlp007 = 'Rev007E_RSS_2005_123_X43_E/RSS_2005_123_X43_E_DLP_500M.TAB'
    csv007 = ExtractCSVData(geo007, cal007, dlp007)
    wtypes = ["KB25", "KBMD25", "KB20", "KBMD20"]
    wp = []
    wf = []
    ym = []
    data = DiffractionCorrection(csv007, 1, 'kb25', bfac = True)
    r = data.rho_km_vals
    xmin = numpy.min(r)
    xmax = numpy.max(r)
    xlim = xmin, xmax
    nc = numpy.min((r >= 100000).nonzero())

    for wtype in wtypes:
        for bfac in [True, False]:
            data = DiffractionCorrection(
                csv007, 1, wtype = wtype, bfac = bfac
            )

            wfnc = data.w_km_vals
            wval = wfnc[nc]
            ymax = numpy.max(wfnc)
            wp.append(wfnc)
            wf.append(wval)
            ym.append(ymax)

    return r, wtypes, wp, wf, ym, xlim

def main():
    """
        Make the plots.
    """
    _, wtypes, _, wf, _, _ = rev007_window()

    date = date_string()
    nwcoms = numpy.size(wf)
    savestring = "rev007_" + date + ".pdf"

    s = []
    s.append(r"\begin{table}")
    s.append(r"\centering")
    s.append(r"\begin{tabular}{|l|l|l|l|}")
    s.append(r"\hline")
    s.append(r"Window Type&Resolution (km)&Allen Deviation&Window Width (km)\\")
    s.append(r"\hline")

    for iwin in range(int(nwcoms/2.0)):
        s.append(r"%s&1.000&$2\times 10^{13})$&%f\\" % (wtypes[iwin], wf[2*iwin]))
        s.append(r"%s&1.000&0.000&%f\\" % (wtypes[iwin], wf[2*iwin + 1]))

    del iwin

    s.append(r"\hline")
    s.append(r"\end{tabular}")
    s.append(r"\caption{Window Widths for Rev007 E X43 at 100,000 km}")
    s.append(r"\end{table}")

    with PdfPages(savestring) as pdf:
        ax = plt.figure(figsize = (8.5, 11))
        ax.text(0.5, 0.5, "outfile.pdf", ha = "center", va = "top")
        pdf.savefig()

    latexstr = ""
    for line in s:
        latexstr = "%s%s" % (latexstr, line)

    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')
    plt.rc('font', size = 10)

    with PdfPages(savestring) as pdf:
        ax = plt.figure(figsize = (8.5, 11))
        ax.text(0.5, 0.5, latexstr, ha = "center", va = "top")
        pdf.savefig()

main()
