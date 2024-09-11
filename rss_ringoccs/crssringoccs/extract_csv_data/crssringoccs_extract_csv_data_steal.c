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

/*  libtmpl provides Booleans and tmpl_String_Duplicate.                      */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>

/*  rssringoccs_CSVData data type typedef'd here.                             */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Function prototype and typedefs for structs given here.                   */
#include "crssringoccs_extract_csv_data.h"

/*  The numpy-array constructor is declared here.                             */
#include "../crssringoccs.h"

/*  Macro for the crssringoccs_set_var function to shorten the syntax.        */
#define CREATE_NUMPY_ARRAY(a) \
crssringoccs_Create_Real_Numpy_Array(&py_csv->a, csv->a, csv->n_elements)

/*  Steals the references to the data in a rssringoccs_CSVData object and     *
 *  creates numpy arrays from them. The data is stored in an instance of the  *
 *  ExtractCSVData class.                                                     */
void
crssringoccs_ExtractCSVData_Steal(crssringoccs_PyCSVObj *py_csv,
                                  rssringoccs_CSVData *csv)
{
    if (!csv || !py_csv)
        return;

    if (csv->error_occurred)
        return;

    if (!py_csv)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message = tmpl_String_Duplicate(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_ExtractCSVData_Steal\n\n"
            "\rInput py_csv is NULL. Aborting.n"
        );

        return;
    }

    CREATE_NUMPY_ARRAY(rho_km_vals);
    CREATE_NUMPY_ARRAY(B_deg_vals);
    CREATE_NUMPY_ARRAY(D_km_vals);
    CREATE_NUMPY_ARRAY(f_sky_hz_vals);
    CREATE_NUMPY_ARRAY(p_norm_vals);
    CREATE_NUMPY_ARRAY(raw_tau_vals);
    CREATE_NUMPY_ARRAY(phase_deg_vals);
    CREATE_NUMPY_ARRAY(phi_deg_vals);
    CREATE_NUMPY_ARRAY(phi_rl_deg_vals);
    CREATE_NUMPY_ARRAY(rho_dot_kms_vals);
    CREATE_NUMPY_ARRAY(t_oet_spm_vals);
    CREATE_NUMPY_ARRAY(t_ret_spm_vals);
    CREATE_NUMPY_ARRAY(t_set_spm_vals);
    CREATE_NUMPY_ARRAY(rx_km_vals);
    CREATE_NUMPY_ARRAY(ry_km_vals);
    CREATE_NUMPY_ARRAY(rz_km_vals);
    CREATE_NUMPY_ARRAY(raw_tau_threshold_vals);
    CREATE_NUMPY_ARRAY(rho_corr_pole_km_vals);
    CREATE_NUMPY_ARRAY(rho_corr_timing_km_vals);
    CREATE_NUMPY_ARRAY(tau_power);
    CREATE_NUMPY_ARRAY(tau_phase);
    CREATE_NUMPY_ARRAY(tau_vals);
}
/*  End of crssringoccs_C_CSV_to_Py_CSV.                                      */

/*  Undefine everything in case someone wants to #include this file.          */
#undef CREATE_NUMPY_ARRAY
