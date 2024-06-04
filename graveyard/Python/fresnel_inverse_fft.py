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
#   Purpose:                                                                   #
#       Computes the numerical Fresnel inverse transform.                      #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/04/17                                                         #
################################################################################
"""
import numpy

def fresnel_inverse_fft(t_hat_vals, kernel, displacement, fresnel_scale):
    """
        Function:
            fresnel_inverse_fft
        Purpose:
            Compute the discrete Fresnel inverse transform via FFT.
        Arguments:
            t_hat_vals:
                The input function (diffraction data).
            kernel:
                The Fresnel kernel.
            displacement:
                The size of a bin (r[1] - r[0]).
            fresnel_scale:
                The Frensel scale.
        Output:
            t_vals:
                The reconstructed data.
        History:
            2024/06/04: Ryan Maguire
                Cleaned up and added to the graveyard directory.
    """
    factor = (0.5 + 0.5j) * displacement / fresnel_scale
    fft_t_hat = numpy.fft.fft(t_hat_vals)
    fft_conv = numpy.fft.fft(kernel)
    inv_t_hat = numpy.fft.ifftshift(numpy.fft.ifft(fft_t_hat*fft_conv))*factor
    return inv_t_hat[(numpy.size(t_hat_vals) - 1) >> 1]
