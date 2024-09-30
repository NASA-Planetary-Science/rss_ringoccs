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
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  The init function for the GetMergedCSVData class. This is the             *
 *  equivalent of the __init__ method defined in a normal python class.       */
int
crssringoccs_GetMergedCSVData_Init(crssringoccs_PyCSVObj *self,
                                PyObject *args,
                                PyObject *kwds)
{
    /*  Variable for the output CSV object.                                   */
    rssringoccs_MergedCSVData *csv = NULL;

    /*  The list of the keywords accepted by the ExtractCSVData class. The    *
     *  file paths 'geo', 'cal', and 'dlp' are required. The 'tau' path is    *
     *  optional, and by default use_deprecate is set to False.               */
    static char *kwlist[] = {"dlpm", NULL};

    /*  Python objects needed throughout the computation.                     */
    PyObject *tmp;

    /*  Initialize the strings for the file paths to NULL. This helps the     *
     *  inner C routines detect errors.                                       */
    const char *dlpm_str = NULL;

    /*  Extract the inputs and keywords supplied by the user. If the data     *
     *  cannot be extracted, raise a type error and return to caller. A short *
     *  explanation of PyArg_ParseTupleAndKeywords. The inputs args and kwds  *
     *  are somewhat straight-forward, they're the arguments and keywords     *
     *  passed by the string. The cryptic string is not straight-forward. The *
     *  | symbol means everything after need not be positional, and we can    *
     *  specify arguments and keywords by name when calling                   *
     *  ExtractCSVData, for example ExtractCSVData(..., tau = "TAU.TAB").     *
     *  s indicates a string, which are the paths to the GEO, CAL, and DLP    *
     *  files. The $ symbol means everything after is optional. p is a        *
     *  Boolean (p for "predicate"), this is the use_deprecate keyword. The   *
     *  colon : denotes that the input list has ended.                        */
    int success = PyArg_ParseTupleAndKeywords(
        args,
        kwds,
        "|s:",
        kwlist,
        &dlpm_str
    );

    if (!success)
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tGetMergedCSVData\n\n"
            "\rCould not parse input variables.\n\n"
            "\rInputs:\n"
            "\r\tdlpm:          Location of a DLPM.TAB file (str)\n"
        );
        return -1;
    }

    /*  All of the heavy lifting is done via the C routines.                  */
    csv = rssringoccs_MergedCSVData_Extract(dlpm_str);

    /*  Several things can fail in this process. Firstly, malloc can fail     *
     *  and return a NULL pointer. Check for this.                            */
    if (!csv)
    {
        PyErr_Format(
            PyExc_RuntimeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tGetMergedCSVData\n\n"
            "\rrssringoccs_CSVData_Extract returned NULL. Aborting.\n"
        );

        return -1;
    }

    /*  Many errors can occur while parsing the data. Check the flag.         */
    if (csv->error_occurred)
    {
        /*  Malloc may have failed to allocate memory for the error message.  *
         *  If so, print a generic error.                                     */
        if (!csv->error_message)
        {
            PyErr_Format(
                PyExc_RuntimeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tGetMergedCSVData\n\n"
                "\rrssringoccs_CSVData_Extract returned with error_occurred\n"
                "\rset to True. No error message was set. Aborting.\n"
            );
        }

        /*  Otherwise print the error message set in the C routine. This can  *
         *  greatly help with debugging.                                      */
        else
        {
            PyErr_Format(
                PyExc_RuntimeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tGetMergedCSVData\n\n"
                "\rrssringoccs_CSVData_Extract returned with error_occurred\n"
                "\rset to True. The following error message was set:\n\n"
                "%s",
                csv->error_message
            );
        }

        /*  Since the function did not return NULL, it is likely that memmory *
         *  was allocated to some variables inside the struct. Free them.     */
        rssringoccs_MergedCSVData_Destroy(&csv);
        return -1;
    }

    /*  To avoid duplicating memory, the Python object simply steals the data *
     *  inside the C object.                                                  */
    crssringoccs_GetMergedCSVData_Steal(self, csv);

    /*  Log how this object was created. Add the history object.              */
    crssringoccs_GetMergedCSVData_Create_History(self, dlpm_str);

    tmp = self->rev_info;
    Py_INCREF(Py_None);
    self->rev_info = Py_None;
    Py_XDECREF(tmp);

    return 1;
}

#undef MAKE_NONE
