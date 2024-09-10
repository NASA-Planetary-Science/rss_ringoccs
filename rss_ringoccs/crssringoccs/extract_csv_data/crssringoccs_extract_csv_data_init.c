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

/*  Function prototype and typedefs for structs given here.                   */
#include "crssringoccs_extract_csv_data.h"

/*  Macro for safely creating None objects.                                   */
#define MAKE_NONE(var)                                                         \
    do {                                                                       \
        PyObject *tmp = self->var;                                             \
        Py_INCREF(Py_None);                                                    \
        self->var = Py_None;                                                   \
        Py_XDECREF(tmp);                                                       \
    } while(0)

/*  The init function for the dirrection correction class. This is the        *
 *  equivalent of the __init__ function defined in a normal python class.     */
int
crssringoccs_ExtractCSVData_Init(crssringoccs_PyCSVObj *self,
                                 PyObject *args,
                                 PyObject *kwds)
{
    rssringoccs_CSVData *csv = NULL;
    tmpl_Bool dpr = tmpl_False;

    /*  The list of the keywords accepted by the DiffractionCorrection class. *
     *  dlp and res are REQUIRED inputs, the rest are optional. If the user   *
     *  does not provide these optional keywords, we must set them ourselves. */
    static char *kwlist[] = {"geo", "cal", "dlp", "use_deprecate", "tau", NULL};

    /*  Python objects needed throughout the computation.                     */
    PyObject *tmp, *csv_tmp;
    const char *geo_str = NULL;
    const char *cal_str = NULL;
    const char *dlp_str = NULL;
    const char *tau_str = NULL;

    /*  Extract the inputs and keywords supplied by the user. If the data     *
     *  cannot be extracted, raise a type error and return to caller. A short *
     *  explaination of PyArg_ParseTupleAndKeywords. The inputs args and kwds *
     *  are somewhat straight-forward, they're the arguments and keywords     *
     *  passed by the string. The cryptic string is not straight-forward. The *
     *  | symbol means everything after need not be positional, and we can    *
     *  specify arguments and keywords by name when calling                   *
     *  DiffractionCorrection, for example                                    *
     *  DiffractionCorrect(..., wtype="blah"). O indicates a Python object,   *
     *  and d is a Python float. This is the DLP and res variables. The $     *
     *  symbold means everything after is optional. s is a string, p is a     *
     *  Boolean (p for "predicate"). b is an integer, and the colon : denotes *
     *  that the input list has ended.                                        */
    int success = PyArg_ParseTupleAndKeywords(
        args,
        kwds,
        "|sss$ps:",
        kwlist,
        &geo_str,
        &cal_str,
        &dlp_str,
        &dpr,
        &tau_str
    );

    if (!success)
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tExtractCSVData\n\n"
            "\rCould not parse input variables.\n\n"
            "\rInputs:\n"
            "\r\tgeo:           Location of a GEO.TAB file (str)\n"
            "\r\tcal:           Location of a CAL.TAB file (str)\n"
            "\r\tdlp:           Location of a DLP.TAB file (str)\n"
            "\rKeywords:\n"
            "\r\tuse_deprecate: Use the old CSV format.\n"
            "\r\ttau:           Location of a TAU.TAB file (str)\n"
        );
        return -1;
    }

    csv = rssringoccs_CSVData_Extract(geo_str, cal_str, dlp_str, tau_str, dpr);

    if (!csv)
    {
        PyErr_Format(
            PyExc_RuntimeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tExtractCSVData\n\n"
            "\rrssringoccs_CSVData_Extract returned NULL. Aborting.\n"
        );

        return -1;
    }

    if (csv->error_occurred)
    {
        if (!csv->error_message)
        {
            PyErr_Format(
                PyExc_RuntimeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tExtractCSVData\n\n"
                "\rrssringoccs_CSVData_Extract returned with error_occurred\n"
                "\rset to True. No error message was set. Aborting.\n"
            );
        }
        else
        {
            PyErr_Format(
                PyExc_RuntimeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tExtractCSVData\n\n"
                "\rrssringoccs_CSVData_Extract returned with error_occurred\n"
                "\rset to True. The following error message was set:\n\n"
                "%s",
                csv->error_message
            );
        }

        rssringoccs_CSVData_Destroy(&csv);
        return -1;
    }

    crssringoccs_ExtractCSVData_Steal(self, csv);

    csv_tmp = Py_BuildValue(
        "{s:s,s:s,s:s}",
        "geo", geo_str,
        "cal", cal_str,
        "dlp", dlp_str
    );

    tmp = self->input_vars;
    Py_INCREF(csv_tmp);
    self->input_vars = csv_tmp;
    Py_XDECREF(tmp);

    tmp = self->input_kwds;
    Py_INCREF(Py_None);
    self->input_kwds = Py_None;
    Py_XDECREF(tmp);

    tmp = self->rev_info;
    Py_INCREF(Py_None);
    self->rev_info = Py_None;
    Py_XDECREF(tmp);

    /*  TODO: Handle history correctly. It's none for now to avoid segfaults. */
    MAKE_NONE(history);

    return 1;
}

#undef MAKE_NONE
