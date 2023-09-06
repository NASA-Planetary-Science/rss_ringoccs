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
 *      Typedef for the Tau object and functions working with the struct.     *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       July 21, 2023                                                 *
 ******************************************************************************/

/*  Include guard to prevent including this file twice.                       */
#ifndef RSS_RINGOCCS_TAU_H
#define RSS_RINGOCCS_TAU_H

/*  Booleans given here.                                                      */
#include <libtmpl/include/tmpl_bool.h>

/*  Complex numbers provided here.                                            */
#include <libtmpl/include/tmpl_complex.h>

/*  size_t typedef provided here.                                             */
#include <stddef.h>

/*  Window function, input is x-parameter and window width.                   */
typedef double (*rssringoccs_Window_Function)(double, double);

/*  All of the types of reconstruction.                                       */
typedef enum {

    /*  Fresnel quadratic approximation.                                      */
    rssringoccs_DR_Fresnel = 0,

    /*  Legendre polynomials, higher order Fresnel approximations.            */
    rssringoccs_DR_Legendre = 1,

    /*  Newton-Raphson method, and various interpolating methods.             */
    rssringoccs_DR_Newton = 2,
    rssringoccs_DR_NewtonQuadratic = 3,
    rssringoccs_DR_NewtonQuartic = 4,
    rssringoccs_DR_NewtonSextic = 5,
    rssringoccs_DR_NewtonOctic = 6,

    /*  Newton-Raphson method with D correction, and interpolations.          */
    rssringoccs_DR_NewtonD = 7,
    rssringoccs_DR_NewtonDQuadratic = 8,
    rssringoccs_DR_NewtonDQuartic = 9,
    rssringoccs_DR_NewtonDSextic = 10,
    rssringoccs_DR_NewtonDOctic = 11,

    /*  Newton-Raphson method with the old D correction, and interpolations.  */
    rssringoccs_DR_NewtonDOld = 12,
    rssringoccs_DR_NewtonDOldQuadratic = 13,
    rssringoccs_DR_NewtonDOldQuartic = 14,
    rssringoccs_DR_NewtonDOldSextic = 15,
    rssringoccs_DR_NewtonDOldOctic = 16,

    /*  Newton-Raphson with dD / dphi correction, and interpolations.         */
    rssringoccs_DR_NewtonDPhi = 17,
    rssringoccs_DR_NewtonDPhiQuadratic = 18,
    rssringoccs_DR_NewtonDPhiQuartic = 19,
    rssringoccs_DR_NewtonDPhiSextic = 20,
    rssringoccs_DR_NewtonDPhiOctic = 21,

    /*  Newton-Raphson with arbitrary quartic perturbation polynomial.        */
    rssringoccs_DR_NewtonPerturb = 22,

    /*  Newton-Raphson, using a single FFT across the entire data set.        */
    rssringoccs_DR_NewtonSimpleFFT = 23,

    /*  Newton-Raphson with elliptical corrections, and interpolations.       */
    rssringoccs_DR_NewtonElliptical = 24,
    rssringoccs_DR_NewtonEllipticalQuadratic = 25,
    rssringoccs_DR_NewtonEllipticalQuartic = 26,
    rssringoccs_DR_NewtonEllipticalSextic = 27,
    rssringoccs_DR_NewtonEllipticalOctic = 28
} rssringoccs_Psitype_Enum;

/*  Structure that contains all of the necessary data.                        */
typedef struct rssringoccs_TAUObj_Def {
    tmpl_ComplexDouble *T_in;
    tmpl_ComplexDouble *T_out;
    tmpl_ComplexDouble *T_fwd;
    double *rho_km_vals;
    double *F_km_vals;
    double *phi_deg_vals;
    double *k_vals;
    double *rho_dot_kms_vals;
    double *B_deg_vals;
    double *D_km_vals;
    double *w_km_vals;
    double *t_oet_spm_vals;
    double *t_ret_spm_vals;
    double *t_set_spm_vals;
    double *rho_corr_pole_km_vals;
    double *rho_corr_timing_km_vals;
    double *tau_threshold_vals;
    double *phi_rl_deg_vals;
    double *rx_km_vals;
    double *ry_km_vals;
    double *rz_km_vals;
    double dx_km;
    double normeq;
    double sigma;
    double ecc;
    double peri;
    double res;
    double perturb[5];
    double rng_list[2];
    double rng_req[2];
    double EPS;
    unsigned int toler;
    size_t start;
    size_t n_used;
    size_t arr_size;
    rssringoccs_Window_Function window_func;
    rssringoccs_Psitype_Enum psinum;
    tmpl_Bool use_norm;
    tmpl_Bool use_fwd;
    tmpl_Bool bfac;
    tmpl_Bool verbose;
    tmpl_Bool error_occurred;
    char *error_message;
    unsigned int order;
} rssringoccs_TAUObj;

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Destroy                                               *
 *  Purpose:                                                                  *
 *      Frees all data associated with a Tau object.                          *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj **):                                          *
 *          The Tau object that is to be destroyed.                           *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_Tau_Destroy(rssringoccs_TAUObj **tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Destroy_Members                                       *
 *  Purpose:                                                                  *
 *      Frees all data associated with a Tau object without destroying the    *
 *      actual Tau pointer.                                                   *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object whose members are to be freed.                     *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_Tau_Destroy_Members(rssringoccs_TAUObj *tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Set_Psi_Type                                          *
 *  Purpose:                                                                  *
 *      Sets the reconstruction method corresponding to a string.             *
 *  Arguments:                                                                *
 *      psitype (const char *):                                               *
 *          A string containing the requested reconstruction method.          *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object whose psinum is to be set.                         *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void
rssringoccs_Tau_Set_Psi_Type(const char *psitype, rssringoccs_TAUObj* tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Set_Window_Type                                       *
 *  Purpose:                                                                  *
 *      Sets the window function in a Tau object corresponding to a string.   *
 *  Arguments:                                                                *
 *      wtype (const char *):                                                 *
 *          A string containing the requested window funtion.                 *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object whose window function is to be set.                *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void
rssringoccs_Tau_Set_Window_Type(const char *wtype, rssringoccs_TAUObj *tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Set_Range_From_String                                 *
 *  Purpose:                                                                  *
 *      Sets the reconstruction range corresponding to a string.              *
 *  Arguments:                                                                *
 *      range (const char *):                                                 *
 *          A string containing the requested range of reconstruction.        *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object whose range is to be set.                          *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void
rssringoccs_Tau_Set_Range_From_String(const char *range,
                                      rssringoccs_TAUObj* tau);

extern void
rssringoccs_Copy_DLP_Data_To_Tau(const rssringoccs_DLPObj *dlp,
                                 rssringoccs_TAUObj *tau);

#endif
/*  End of include guard.                                                     */
