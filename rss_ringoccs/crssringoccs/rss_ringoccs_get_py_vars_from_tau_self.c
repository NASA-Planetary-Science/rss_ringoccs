/******************************************************************************
 *                                 LICENSE                                    *
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
#include <rss_ringoccs/include/crss_ringoccs.h>

void
rssringoccs_Get_Py_Vars_From_Tau_Self(rssringoccs_TAUObj *tau,
                                      PyDiffrecObj *self)
{
    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    if (self == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_C_Tau_to_Py_Tau\n\n"
            "\rInput self is NULL. Aborting.n"
        );
        return;
    }

    tau->sigma    = self->sigma;
    tau->bfac     = self->bfac;
    tau->ecc      = self->ecc;
    tau->peri     = self->peri;
    tau->use_fwd  = self->use_fwd;
    tau->use_norm = self->use_norm;
    tau->verbose  = self->verbose;
}

