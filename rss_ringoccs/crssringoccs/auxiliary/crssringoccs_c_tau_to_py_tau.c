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

/*  NULL is defined here.                                                     */
#include <stddef.h>

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  tmpl_strdup function declared here.                                       */
#include <libtmpl/include/tmpl_string.h>

/*  Function prototype and typedefs for structs given here.                   */
#include "../crssringoccs.h"

/*  Macro for the crssringoccs_set_var function to shorten the syntax.        */
#define SET_VAR(a) crssringoccs_Set_Var(&py_tau->a, tau->a, tau->arr_size)
#define SET_CVAR(a) crssringoccs_Set_CVar(&py_tau->a, tau->a, tau->arr_size)

/*  Macro for safely creating None objects.                                   */
#define MAKE_NONE(var)                                                         \
    do {                                                                       \
        PyObject *tmp = py_tau->var;                                           \
        Py_INCREF(Py_None);                                                    \
        py_tau->var = Py_None;                                                 \
        Py_XDECREF(tmp);                                                       \
    } while(0)

/*  Converts a C Tau struct to a Python Tau Object.                           */
void crssringoccs_C_Tau_To_Py_Tau(PyDiffrecObj *py_tau, rssringoccs_TAUObj *tau)
{
    /*  If the C version of the object is NULL there is nothing to do.        */
    if (tau == NULL)
        return;

    /*  Do not attempt to convert if an error occurred while making tau.      */
    if (tau->error_occurred)
        return;

    /*  If the Python object hasn't been allocated memory, abort with error.  */
    if (py_tau == NULL)
    {
        tau->error_occurred = tmpl_True;

        /*  Create an error message so the user knows what went wrong.        */
        tau->error_message = tmpl_String_Duplicate(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_C_Tau_To_Py_Tau\n\n"
            "\rInput py_tau is NULL. Aborting.n"
        );
        return;
    }

    /*  Set every variable in the Python object from the C Tau struct.        */
    SET_CVAR(T_in);
    SET_CVAR(T_out);
    SET_VAR(rho_km_vals);
    SET_VAR(F_km_vals);
    SET_VAR(phi_deg_vals);
    SET_VAR(k_vals);
    SET_VAR(rho_dot_kms_vals);
    SET_VAR(B_deg_vals);
    SET_VAR(D_km_vals);
    SET_VAR(w_km_vals);
    SET_VAR(t_oet_spm_vals);
    SET_VAR(t_ret_spm_vals);
    SET_VAR(t_set_spm_vals);
    SET_VAR(rho_corr_pole_km_vals);
    SET_VAR(rho_corr_timing_km_vals);
    SET_VAR(tau_threshold_vals);
    SET_VAR(phi_rl_deg_vals);
    SET_VAR(rx_km_vals);
    SET_VAR(ry_km_vals);
    SET_VAR(rz_km_vals);

    /*  If forward modeling was not performed, set these as None objects.     */
    if (tau->T_fwd == NULL)
        MAKE_NONE(T_fwd);

    else
        SET_CVAR(T_fwd);
}
/*  End of crssringoccs_C_Tau_To_Py_Tau.                                      */

/*  Undefine these in case someone wants to #include this file.               */
#undef SET_VAR
#undef MAKE_NONE
