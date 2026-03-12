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
#include "../crssringoccs.h"

void
crssringoccs_DiffractionCorrection_Create_Argument_Dictionary(
    crssringoccs_PyDiffrecObj * const self,
    PyObject * const dlp
)
{
    PyObject *tmp = NULL;
    PyObject *input_vars = NULL;
    PyObject *dlp_history = NULL;

    if (!self)
        return;

    if (!self->tau)
        return;

    if (self->tau->error_occurred)
        return;

    if (self->verbose)
        puts("\r\tDiffractionCorrection: Creating argument dictionary...");

    if (!dlp)
    {
        self->tau->error_occurred = tmpl_True;
        self->tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_DiffractionCorrection_Create_Argument_Dictionary"
            "\n\n\rInput DLP object is NULL.\n\n";

        return;
    }

    if (!PyObject_HasAttrString(dlp, "history"))
    {
        self->tau->error_occurred = tmpl_True;
        self->tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_DiffractionCorrection_Create_Argument_Dictionary"
            "\n\n\rInput DLP object does not have a 'history' attribute.\n\n";

        return;
    }

    dlp_history = PyObject_GetAttrString(dlp, "history");

    if (!dlp_history)
    {
        self->tau->error_occurred = tmpl_True;
        self->tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_DiffractionCorrection_Create_Argument_Dictionary"
            "\n\n\rPyObject_GetAttrString returned NULL for 'history'.\n\n";

        return;
    }

    input_vars = Py_BuildValue(
        /*  Format specifier for the argument dictionary. This produces:      *
         *      String: Python object                                         *
         *      String: double                                                *
         *  The arguments are the DLP object and the resolution.              */
        "{s:O,s:d}",

        /*  The two members of the dictionary, the DLP and the resolution.    */
        "dlp_inst", dlp_history,
        "resolution_km", self->input_resolution_km
    );

    if (!input_vars)
    {
        self->tau->error_occurred = tmpl_True;
        self->tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_DiffractionCorrection_Create_Argument_Dictionary"
            "\n\n\rPy_BuildValue returned NULL.\n\n";

        Py_CLEAR(dlp_history);
        return;
    }

    tmp = self->input_vars;
    self->input_vars = input_vars;

    Py_CLEAR(tmp);
    Py_CLEAR(dlp_history);
}
/*  End of crssringoccs_DiffractionCorrection_Create_Argument_Dictionary.     */
