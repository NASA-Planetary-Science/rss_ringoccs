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

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  Function prototype and typedefs for structs given here.                   */
#include "../crssringoccs.h"

void
crssringoccs_Get_Py_Vars_From_Tau_Self(rssringoccs_TAUObj *tau,
                                       const crssringoccs_PyDiffrecObj *self)
{
    if (!tau)
        return;

    if (tau->error_occurred)
        return;

    if (!self)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_Get_Py_Vars_From_Tau_Self\n\n"
            "\rInput self is NULL.n";

        return;
    }

    tau->sigma = self->sigma;
    tau->bfac = self->bfac;
    tau->eccentricity = self->eccentricity;
    tau->periapse = self->periapse;
    tau->use_fwd = self->use_fwd;
    tau->use_norm = self->use_norm;
    tau->verbose = self->verbose;
}
