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
#include "../crssringoccs.h"

void
crssringoccs_DiffractionCorrection_Set_Range(
    crssringoccs_PyDiffrecObj * const self
)
{
    PyObject *iter;
    PyObject *next;
    unsigned int n;

    if (!self)

    if (!self->tau)
        return;

    if (self->tau->error_occurred)
        return;

    if (self->verbose)
        puts("\r\tDiffractionCorrection: Passing 'rng' to C struct...");

    if (!self->rngreq)
    {
        self->rngreq = PyUnicode_FromString("all");

        if (!self->rngreq)
        {
            self->tau->error_occurred = tmpl_True;
            self->tau->error_message =
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tcrssringoccs_DiffractionCorrection_Set_Range\n\n"
                "\rInput rngreq is NULL. Aborting.\n\n";

            return;
        }
    }

    /*  If the rng variable is a string, make sure it is a legal value and    *
     *  try to extract the corresponding values in kilometers.                */
    if (PyBytes_Check(self->rngreq))
    {
        const char * const range_string = PyBytes_AsString(self->rngreq);
        rssringoccs_Tau_Set_Range_From_String(range_string, self->tau);
    }

    /*  If the rng variable is a unicode object (type of string from python)  *
     *  make sure it is a legal value and try to extract the corresponding    *
     *  values in kilometers.                                                 */
    else if (PyUnicode_Check(self->rngreq))
    {

        /*  Convert the Python string to a C string via PyUnicode_AsUTF8. The *
         *  C API recommends not altering the string, so we create a copy of  *
         *  it using strcpy (from string.h).                                  */
        const char * const range_string = PyUnicode_AsUTF8(self->rngreq);
        rssringoccs_Tau_Set_Range_From_String(range_string, self->tau);
    }

    /*  If the requested range is a list, try to parse the elements.          */
    else if (PyList_Check(self->rngreq))
    {
        if (PyList_Size(self->rngreq) != 2)
        {
            self->tau->error_occurred = tmpl_True;
            self->tau->error_message =
                "\rError Encountered: rss_ringoccs\n"
                "\r\tcrssringoccs_DiffractionCorrection_Set_Range\n\n"
                "\rInput range is a list but does not have 2 entries.\n"
                "\rrng must be a list of two real numbers.\n\n";

            return;
        }

        iter = PyObject_GetIter(self->rngreq);

        for (n = 0; n < 2; ++n)
        {
            next = PyIter_Next(iter);

            /*  Try to parse the elements. Return with error if this fails.   */
            if (PyLong_Check(next))
                self->tau->range[n] = PyLong_AsDouble(next);
            else if (PyFloat_Check(next))
                self->tau->range[n] = PyFloat_AsDouble(next);
            else
            {
                self->tau->error_occurred = tmpl_True;
                self->tau->error_message =
                    "\rError Encountered: rss_ringoccs\n"
                    "\r\tcrssringoccs_DiffractionCorrection_Set_Range\n\n"
                    "\rInput rng has entries that are not real numbers.\n"
                    "\rBoth entries for the rng list must be numbers.\n\n";

                return;
            }

            Py_CLEAR(next);
        }

        Py_CLEAR(iter);
    }

    /*  Illegal rng requested. Return with error.                             */
    else
    {
        self->tau->error_occurred = tmpl_True;
        self->tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Get_Py_Range\n\n"
            "\rrng must be a list of two real numbers or a string.\n"
            "\rAllowed strings are:\n"
            "\r\tall               [1.0, 400000.0]\n"
            "\r\tbesselbarnard     [120210.0, 120330.0]\n"
            "\r\tbessel-barnard    [120210.0, 120330.0]\n"
            "\r\tcringripples      [77690.0, 77760.0]\n"
            "\r\tencke             [132900.0, 134200.0]\n"
            "\r\tenckegap          [132900.0, 134200.0]\n"
            "\r\therschel          [118100.0, 118380.0]\n"
            "\r\therschelgap       [118100.0, 118380.0]\n"
            "\r\thuygens           [117650.0, 117950.0]\n"
            "\r\thuygensringlet    [117650.0, 117950.0]\n"
            "\r\tjanusepimetheus   [96200.0, 96800.0]\n"
            "\r\tjeffreys          [118900.0, 119000.0]\n"
            "\r\tjeffreysgap       [118900.0, 119000.0]\n"
            "\r\tkuiper            [119300.0, 119500.0]\n"
            "\r\tkuipergap         [119300.0, 119500.0]\n"
            "\r\tmaxwell           [87410.0, 87610.0]\n"
            "\r\tmaxwellringlet    [87410.0, 87610.0]\n"
            "\r\trussell           [118550.0, 118660.0]\n"
            "\r\trussellgap        [118550.0, 118660.0]\n"
            "\r\ttitan             [77870.0, 77930.0]\n"
            "\r\ttitanringlet      [77870.0, 77930.0]\n\n";

        return;
    }
}
