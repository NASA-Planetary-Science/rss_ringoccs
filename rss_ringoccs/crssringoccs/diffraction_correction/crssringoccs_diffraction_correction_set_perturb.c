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
crssringoccs_DiffractionCorrection_Set_Perturb(
    crssringoccs_PyDiffrecObj * const self
)
{
    PyObject *iter;
    PyObject *next;
    unsigned int n;

    if (!self)
        return;

    if (!self->tau)
        return;

    if (self->tau->error_occurred)
        return;

    if (self->verbose)
        puts("\r\tDiffractionCorrection: Passing 'perturb' to C struct...");

    /*  Check that the input perturb is a list with 5 elements.               */
    if (!self->perturb)
    {
        for (n = 0; n < 5; ++n)
            self->tau->perturb[n] = 0.0;
    }

    /*  If the user supplied a perturb list, parse it and extract values.     */
    else if (PyList_Check(self->perturb))
    {
        /*  If the list is not the correct size, raise an error.              */
        if (PyList_Size(self->perturb) != 5)
        {
            self->tau->error_occurred = tmpl_True;
            self->tau->error_message =
                "\rError Encountered: rss_ringoccs\n"
                "\r\tcrssringoccs_DiffractionCorrection_Set_Perturb\n\n"
                "\rInput perturb is a list but does not have 5 entries.\n"
                "\rperturb must be a list of five real numbers.\n";

            return;
        }

        iter = PyObject_GetIter(self->perturb);

        /*  Loop over the elements of the list, see if they can be converted  *
         *  to doubles, and store them in the tau->perturb variable.          */
        for (n = 0; n < 5; ++n)
        {
            next = PyIter_Next(iter);

            /*  If the element is an integer, convert to double and save it.  */
            if (PyLong_Check(next))
                self->tau->perturb[n] = PyLong_AsDouble(next);

            /*  Convert from Python float to C double with PyFloat_AsDouble.  */
            else if (PyFloat_Check(next))
                self->tau->perturb[n] = PyFloat_AsDouble(next);

            /*  Invalid data type for one of the entries. Return with error.  */
            else
            {
                self->tau->error_occurred = tmpl_True;
                self->tau->error_message =
                    "\rError Encountered: rss_ringoccs\n"
                    "\r\tcrssringoccs_DiffractionCorrection_Set_Perturb\n\n"
                    "\rInput perturb has entries that are not real numbers.\n"
                    "\rAll entries for the perturb list must be numbers.\n";

                return;
            }

            Py_CLEAR(next);
        }

        Py_CLEAR(iter);
    }

    /*  The input was not a list. Return with error.                          */
    else
    {
        self->tau->error_occurred = tmpl_True;
        self->tau->error_message =
            "\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_DiffractionCorrection_Set_Perturb\n\n"
            "\rInput perturb is not a list.\n";

        return;
    }
}
