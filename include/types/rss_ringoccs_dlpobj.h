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
 *      Typedef for the DLP object with geometry and diffraction data.        *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       February 27, 2026                                             *
 ******************************************************************************/

/*  Include guard to prevent including this file twice.                       */
#ifndef RSS_RINGOCCS_TYPES_DLPOBJ_H
#define RSS_RINGOCCS_TYPES_DLPOBJ_H

/*  Booleans found here.                                                      */
#include <libtmpl/include/tmpl_bool.h>

/*  size_t typedef given here.                                                */
#include <stddef.h>

/*  Structure that contains all of the necessary diffraction limited data.    *
 *  This includes geometry data, frequency data, and the diffraction profile. */
typedef struct rssringoccs_DLPObj_Type {
    double *rho_km_vals;
    double *phi_deg_vals;
    double *B_deg_vals;
    double *D_km_vals;
    double *f_sky_hz_vals;
    double *rho_dot_kms_vals;
    double *t_oet_spm_vals;
    double *t_ret_spm_vals;
    double *t_set_spm_vals;
    double *rho_corr_pole_km_vals;
    double *rho_corr_timing_km_vals;
    double *phi_rl_deg_vals;
    double *p_norm_vals;
    double *phase_deg_vals;
    double *raw_tau_threshold_vals;
    double *rx_km_vals;
    double *ry_km_vals;
    double *rz_km_vals;
    double dx_km;
    size_t arr_size;

    /*  Other objects include this object as a member. We use reference       *
     *  counting to avoid duplicating the data. This object is destroyed once *
     *  the reference counter hits zero.                                      */
    size_t reference_count;
    tmpl_Bool error_occurred;
    const char *error_message;
} rssringoccs_DLPObj;

#endif
/*  End of include guard.                                                     */
