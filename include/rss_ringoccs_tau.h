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

/*  Booleans and complex numbers given here.                                  */
#include <libtmpl/include/tmpl.h>

/*  size_t typedef provided here.                                             */
#include <stddef.h>

/*  DLP object is typedef'd here.                                             */
#include <rss_ringoccs/include/rss_ringoccs_calibration.h>

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
    rssringoccs_DR_NewtonEllipticalOctic = 28,

    /*  Indicates an error.                                                   */
    rssringoccs_DR_None = 100
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
 *      rssringoccs_Tau_Create_From_DLP                                       *
 *  Purpose:                                                                  *
 *      Creates a Tau object using the data in a DLP objects.                 *
 *  Arguments:                                                                *
 *      dlp (const rssringoccs_DLPObj *):                                     *
 *          The DLP object whose members are being copied.                    *
 *      res (double):                                                         *
 *          The requested resolution, in kilometers.                          *
 *  Outputs:                                                                  *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The output Tau object.                                            *
 *  Notes:                                                                    *
 *      malloc is called to allocate memory for the output pointer and the    *
 *      members in the Tau object. If malloc fails to allocate memory for the *
 *      tau pointer, then NULL is returned. If malloc fails to allocate       *
 *      memory for any of the members of the Tau objects, the error_occurred  *
 *      Boolean is set to true. Check both of these before accessing data.    *
 ******************************************************************************/
extern rssringoccs_TAUObj *
rssringoccs_Tau_Create_From_DLP(const rssringoccs_DLPObj *dlp, double res);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Copy_DLP_Data                                         *
 *  Purpose:                                                                  *
 *      Copies data and computes variables for a Tau object from DLP data.    *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object. DLP members will be copied here.                  *
 *      dlp (const rssringoccs_DLPObj *):                                     *
 *          The DLP object whose members are being copied.                    *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void
rssringoccs_Tau_Copy_DLP_Data(const rssringoccs_DLPObj *dlp,
                              rssringoccs_TAUObj *tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Copy_DLP_Members                                      *
 *  Purpose:                                                                  *
 *      Copies members that dlp and tau objects have in common from a dlp     *
 *      object in to a tau object.                                            *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object. DLP members will be copied here.                  *
 *      dlp (const rssringoccs_DLPObj *):                                     *
 *          The DLP object whose members are being copied.                    *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void
rssringoccs_Tau_Copy_DLP_Members(rssringoccs_TAUObj *tau,
                                 const rssringoccs_DLPObj *dlp);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Compute_Data_From_DLP_Members                         *
 *  Purpose:                                                                  *
 *      Computes several Tau variables from the given DLP data.               *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object.                                                   *
 *      dlp (const rssringoccs_DLPObj *):                                     *
 *          The DLP object.                                                   *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void
rssringoccs_Tau_Compute_Data_From_DLP_Members(rssringoccs_TAUObj *tau,
                                              const rssringoccs_DLPObj *dlp);

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
 *      rssringoccs_Tau_Init                                                  *
 *  Purpose:                                                                  *
 *      Initialize a tau struct so that it's members are NULL.                *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object whose members are to initialized.                  *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 *  Notes:                                                                    *
 *      Functions that allocate memory for a Tau object check if the members  *
 *      are NULL before working with them. If the pointer is not NULL it is   *
 *      assumed memory was successfully allocated for the variable. Always    *
 *      call this function on a new tau object first.                         *
 ******************************************************************************/
extern void rssringoccs_Tau_Init(rssringoccs_TAUObj *tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Malloc_Members                                        *
 *  Purpose:                                                                  *
 *      Allocates memory for the Tau variables.                               *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object whose members are to be allocated memory.          *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 *  Notes:                                                                    *
 *      It is assumed tau->arr_size has been set and that tau is not NULL.    *
 ******************************************************************************/
extern void rssringoccs_Tau_Malloc_Members(rssringoccs_TAUObj *tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Has_Errors                                            *
 *  Purpose:                                                                  *
 *      Checks a Tau object for possible errors.                              *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object to be checked.                                     *
 *  Outputs:                                                                  *
 *      has_errors (tmpl_Bool):                                               *
 *          A Boolean. True if the Tau object contains errors. Otherwise the  *
 *          function returns False.                                           *
 ******************************************************************************/
extern tmpl_Bool rssringoccs_Tau_Has_Errors(rssringoccs_TAUObj *tau);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Tau_Set_Default_Values                                    *
 *  Purpose:                                                                  *
 *      Sets the default values for several parameters.                       *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj *):                                           *
 *          The Tau object whose values are to be set.                        *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_Tau_Set_Default_Values(rssringoccs_TAUObj* tau);

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

#endif
/*  End of include guard.                                                     */
