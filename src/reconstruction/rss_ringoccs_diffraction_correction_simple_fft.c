/******************************************************************************
 *                                  LICENSE                                   *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  rss_ringoccs is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************/
#include <libtmpl/include/compat/tmpl_cast.h>
#include <libtmpl/include/compat/tmpl_free.h>
#include <libtmpl/include/compat/tmpl_malloc.h>
#include <libtmpl/include/tmpl_complex.h>
#include <libtmpl/include/tmpl_math.h>
#include <libtmpl/include/tmpl_fft.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <stddef.h>

static size_t next_pow_2(size_t n)
{
    size_t power = 1;

    while (power < n)
        power <<= 1;

    return power;
}

void rssringoccs_Diffraction_Correction_SimpleFFT(rssringoccs_TAUObj *tau)
{
    /*  Variable for indexing for-loops.                                      */
    size_t n;

    /*  nw_pts is the number of points in the window of integration.          */
    size_t nw_pts, half_nw_pts;

    /*  Variables for indexing the data and processing region.                */
    size_t center, shift, start, data_size, arr_size;

    /*  Some variables needed for reconstruction.                             */
    double w_init, w_thresh;
    tmpl_ComplexDouble *T_in, *ker, *T_out;
    tmpl_ComplexDouble *fft_in, *fft_ker, *fft_out;

    /*  Check that the pointers to the data are not NULL.                     */
    rssringoccs_Tau_Check_Core_Data(tau);

    /* Check to ensure you have enough data around the central data point.    */
    rssringoccs_Tau_Check_Data_Range(tau);

    /*  The previous functions set the error_occurred Boolean on failure.     */
    if (tau->error_occurred)
        return;

    /*  Variable for the center of the data set. The geometry is computed     *
     *  with respect to this point.                                           */
    center = tau->start + tau->n_used / 2;

    /*  Window width needed for the central point, and the number of points   *
     *  in a single window.                                                   */
    w_init = tau->w_km_vals[center];
    half_nw_pts = TMPL_CAST(w_init / tau->dx_km, size_t);
    nw_pts = 2 * half_nw_pts + 1;

    /*  The amount of data that the FFTs process is given by the total number *
     *  of points in the data, plus the half-windows to either side.          */
    data_size = tau->n_used + nw_pts - 1;
    arr_size = next_pow_2(2 * data_size - 1);

    /*  Shift factor for converting the FFT data back into natural order.     */
    shift = data_size >> 1;

    /*  The integral needs half a window to the left of the starting point,   *
     *  and half a window to the right of the end point. Shift the start back.*/
    start = center - shift;

    /*  Allocate memory for the Fresnel kernel and the FFT of the input       *
     *  complex diffracted data.                                              */
    ker = TMPL_MALLOC(tmpl_ComplexDouble, arr_size);
    T_in = TMPL_MALLOC(tmpl_ComplexDouble, arr_size);

    /*  Check if either call to malloc failed.                                */
    if (!ker || !T_in)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Diffraction_Correction_SimpleFFT\n\n"
            "\rmalloc returned NULL.\n";

        /*  ker and fft_in are either NULL, or point to valid data allocated  *
         *  with malloc. In either case we may pass the pointers to free.     */
        TMPL_FREE(ker);
        TMPL_FREE(T_in);
        return;
    }

    /*  We can reuse the memory allocated to the in variables for the out     *
     *  variables. This saves us calls to malloc and wasteful memory use.     */
    fft_ker = ker;
    fft_in = T_in;
    fft_out = fft_ker;
    T_out = fft_out;

    /*  Threshold for computing the Fresnel kernel. Outside of this region    *
     *  the tapering function is zero and kills off the Fresnel kernel, so    *
     *  there is no point in computing it.                                    */
    w_thresh = 0.5 * tau->w_km_vals[center];

    /*  We compute the Fresnel transform using IFFT(FFT(T) * FFT(kernel)).    *
     *  The Tau object already has T, compute the kernel across the data set. *
     *  Note that we are using the geometry for the point in the center of    *
     *  the data set and assuming this approximates the geometry for the      *
     *  entire region being processed. Because of this, trying to process a   *
     *  large region will result in a very poor reconstruction.               */
    for (n = 0; n < data_size; ++n)
    {
        /*  Index for the current point in the window.                        */
        const size_t offset = start + n;

        /*  The Fresnel phase can be expanded in term of rho - rho0. This     *
         *  expression is also the input for the window function.             */
        const double x = tau->rho_km_vals[offset] - tau->rho_km_vals[center];

        /*  Once |rho - rho| exceeds the window width, the output is zero.    */
        if (tmpl_Double_Abs(x) <= w_thresh)
        {
            /*  Compute the product w(x) * exp(i psi(x)) * dx.                */
            const double taper = tau->window_func(x, w_init);
            ker[n] = rssringoccs_Fresnel_Kernel(tau, offset, center);
            tmpl_CDouble_MultiplyBy_Real(&ker[n], taper * tau->dx_km);
        }

        else
            ker[n] = tmpl_CDouble_Zero;

        T_in[n] = tau->T_in[offset];
    }

    for (n = data_size; n < arr_size; ++n)
    {
        ker[n] = tmpl_CDouble_Zero;
        T_in[n] = tmpl_CDouble_Zero;
    }

    /*  By the convolution theorem, FFT(conv(T, ker)) = FFT(T) * FFT(ker).    *
     *  Compute the FFTs of the complex diffraction data and the kernel.      */
    tmpl_CDouble_FFT(T_in, fft_in, arr_size);
    tmpl_CDouble_FFT(ker, fft_ker, arr_size);

    /*  Loop through and compute the point-wise product FFT(T) * FFT(ker).    *
     *  Note, we set fft_out to point to the same data as fft_ker, so this    *
     *  computation reads fft_out = fft_ker * fft_in, or fft_ker *= fft_in.   */
    for (n = 0; n < arr_size; ++n)
        tmpl_CDouble_MultiplyBy(&fft_out[n], &fft_in[n]);

    /*  The original convulation can be obtained by using the inverse FFT on  *
     *  the pointwise product of the two FFTs.                                */
    tmpl_CDouble_IFFT(fft_out, T_out, arr_size);

    /*  Shift the data back so that is in natural order.                      */
    for(n = 0; n < tau->n_used; ++n)
    {
        const size_t i_shift = (half_nw_pts + n + shift) % arr_size;
        tau->T_out[tau->start + n] = T_out[i_shift];
    }

    /*  Free variables allocated by malloc.                                   */
    TMPL_FREE(ker);
    TMPL_FREE(T_in);
}
/*  End of rssringoccs_Diffraction_Correction_SimpleFFT.                      */
