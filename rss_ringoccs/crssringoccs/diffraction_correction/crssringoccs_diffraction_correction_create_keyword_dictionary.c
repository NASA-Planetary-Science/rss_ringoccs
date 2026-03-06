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
crssringoccs_DiffractionCorrection_Create_Keyword_Dictionary(
    crssringoccs_PyDiffrecObj * const self
)
{
    PyObject *input_kwds = NULL;
    PyObject *tmp = NULL;
    PyObject *normalize = NULL;
    PyObject *bfac = NULL;

    if (!self)
        return;

    if (!self->tau)
        return;

    if (self->tau->error_occurred)
        return;

    if (self->verbose)
        puts("\r\tDiffractionCorrection: Creating keyword dictionary...");

    normalize = PyBool_FromLong(self->use_norm);
    bfac = PyBool_FromLong(self->bfac);

    if (!normalize || !bfac)
    {
        self->tau->error_occurred = tmpl_True;
        self->tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_DiffractionCorrection_Create_Keyword_Dictionary"
            "\n\n\rPyBool_FromLong returned NULL.\n";

        Py_CLEAR(normalize);
        Py_CLEAR(bfac);
        return;
    }

    input_kwds = Py_BuildValue(
        "{s:O,s:s,s:s,s:d,s:d,s:d,s:d,s:O,s:O}",
        "rng",               self->rngreq,
        "wtype",             self->wtype,
        "psitype",           self->psitype,
        "sigma",             self->sigma,
        "eccentricity",      self->eccentricity,
        "periapse",          self->periapse,
        "resolution_factor", self->resolution_factor,
        "use_norm",          normalize,
        "bfac",              bfac
    );

    if (!input_kwds)
    {
        self->tau->error_occurred = tmpl_True;
        self->tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_DiffractionCorrection_Create_Keyword_Dictionary"
            "\n\n\rPy_BuildValue returned NULL.\n\n";

        Py_CLEAR(dlp_history);
        return;
    }


    tmp = self->input_kwds;
    Py_INCREF(input_kwds);
    self->input_kwds = input_kwds;

    Py_CLEAR(tmp);
    Py_CLEAR(normalize);
    Py_CLEAR(bfac);
}
/*  End of crssringoccs_DiffractionCorrection_Create_Keyword_Dictionary.      */
