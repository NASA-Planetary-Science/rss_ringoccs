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
#   Author:     Ryan Maguire                                                   #
#   Date:       2018                                                           #
################################################################################
"""
# Pylint wants certain variables to be all caps. The original TC2017
# coding style didn't follow this. Ignore this warning.
# pylint: disable = invalid-name
# pylint: disable = consider-using-f-string
import numpy
from rss_ringoccs import DiffractionCorrection, ExtractCSVData
from .fresnel_kernel import fresnel_kernel

volx = "Archived_Cassini_RSS_RingOccs_2018/data/"
dirx = "Rev133/Rev133E/"
geox = volx + dirx + "RSS_2010_170_X43_E_GEO.TAB"
calx = volx + dirx + "RSS_2010_170_X43_E_CAL.TAB"
dlpx = volx + dirx + "RSS_2010_170_X43_E_DLP_500M.TAB"
taux = volx + dirx + "RSS_2010_170_X43_E_TAU_01KM.TAB"

datax = ExtractCSVData(geox, calx, dlpx, tau=taux, verbose=True)
recx = DiffractionCorrection(
    datax, 0.75, wtype = "kbmd20", verbose = True,
    rng = "cringripples", psitype = "full"
)

# For varying D.
kD = numpy.mean(recx.kD)
r = recx.r
r0 = recx.r0
phi0 = recx.phi0
b = numpy.mean(recx.b)
d = numpy.mean(recx.d)
k = kD/d
psi0 = recx.psi
phi0 = numpy.mean(phi0)
phi = numpy.arange(phi0-0.16, phi0+0.16, 0.00001)
d = numpy.arange(d-1300, d+1300,50.0)
D_arr = numpy.zeros(numpy.size(r0))
dpsi = numpy.zeros(numpy.size(r0))
print("D Calculation")

for i in range(numpy.size(r0)):
    r1 = r0[i]
    psi = fresnel_kernel(
        k*d[None, :], r, r1, phi[:, None], phi0, b, d[None, :]
    )

    psi_min = numpy.amin(psi, axis=0)
    dpsi[i] = numpy.min(numpy.abs(psi_min - psi0[i]))
    D_arr[i] = numpy.min(
        (numpy.abs(psi_min - psi0[i]) == dpsi[i]).nonzero()
    )

    print("\t%d" % (i), end = "\r")

# For varying B.
kD = numpy.mean(recx.kD)
r = recx.r
r0 = recx.r0
phi0 = recx.phi0
b = numpy.mean(recx.b)
d = numpy.mean(recx.d)
k = kD/d
psi0 = recx.psi
phi0 = numpy.mean(phi0)
phi = numpy.arange(phi0-0.16, phi0+0.16, 0.00001)
b = numpy.arange(b-0.0005,b+0.0005,0.00001)
B_arr = numpy.zeros(numpy.size(r0))
dpsi = numpy.zeros(numpy.size(r0))
print("B Calculation")

for i in range(numpy.size(r0)):
    r1 = r0[i]
    psi = fresnel_kernel(kD, r, r1, phi[:, None], phi0, b[None, :], d)
    psi_min = numpy.amin(psi, axis = 0)
    dpsi[i] = numpy.min(numpy.abs(psi_min - psi0[i]))
    B_arr[i] = numpy.min(
        (numpy.abs(psi_min - psi0[i]) == dpsi[i]).nonzero()
    )

    print("\t%d" % (i), end = "\r")

# For varying phi0
kD = numpy.mean(recx.kD)
r = recx.r
r0 = recx.r0
phi0 = recx.phi0
b = numpy.mean(recx.b)
d = numpy.mean(recx.d)
k = kD / d
psi0 = recx.psi
phi0 = numpy.mean(phi0)
phi = numpy.arange(phi0-0.16,phi0+0.16,0.000001)
phi1 = numpy.arange(phi0-0.0005,phi0+0.0005,0.0001)
dpsi = numpy.zeros(numpy.size(r0))
phi_arr = numpy.zeros(numpy.size(r0))
print("phi0 Calculation")

for i in range(numpy.size(r0)):
    r1 = r0[i]
    psi = fresnel_kernel(kD, r, r1, phi[:, None], phi1[None, :], b, d)
    psi_min = numpy.amin(psi, axis=0)
    dpsi[i] = numpy.min(numpy.abs(psi_min - psi0[i]))

    phi_arr[i] = numpy.min((
        numpy.abs(psi_min - psi0[i]) == dpsi[i]).nonzero()
    )

    print("\t%d" % (i), end="\r")
