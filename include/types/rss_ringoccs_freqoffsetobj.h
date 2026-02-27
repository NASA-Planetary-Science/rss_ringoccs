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
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Typedef for the frequency offset object.                              *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       July 21, 2023                                                 *
 ******************************************************************************/

/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_FREQOFFSETOBJ_H
#define RSS_RINGOCCS_FREQOFFSETOBJ_H

/*  Complex numbers found here.                                               */
#include <libtmpl/include/types/tmpl_complex_double.h>

/*  Frequency offset type used for calibration.                               */
typedef struct rssringoccs_FreqOffsetObj_Type {

    /*  Uncorrected real and imaginary components of signal.                  */
    tmpl_ComplexDouble *IQ_m;

    /*  Observed event time at full sampling.                                 */
    double *t_oet_spm_vals;

    /*  Observed event time for frequency.                                    */
    double *t_oet_spm_vals;

    /*  Frequency offset, or frequency at max power.                          */
    double *f_offset_hz_vals;

    /*  Raw time sampling from spm_vals.                                      */
    double dt;

    /*  Half the width of the FFT window.                                     */
    double dt_freq;

    /*  Minimum time for sampling.                                            */
    double t_oet_spm_min;

    /*  Maximum time for sampling.                                            */
    double t_oet_spm_max;
} rssringoccs_FreqOffsetObj;

#endif
/*  End of include guard.                                                     */
