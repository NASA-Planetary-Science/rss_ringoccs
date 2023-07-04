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

void rssringoccs_Get_Py_Perturb(rssringoccs_TAUObj *tau, PyObject *perturb)
{
    PyObject *iter;
    PyObject *next;
    unsigned int n;

    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    /*  Check that the input perturb is a list with 5 elements.               */
    if (!perturb)
    {
        for (n=0; n<5; ++n)
            tau->perturb[n] = 0.0;
    }

    /*  If the user supplied a perturb list, parse it and extract values.     */
    else if (PyList_Check(perturb))
    {
        /*  If the list is not the correct size, raise an error.              */
        if (PyList_Size(perturb) != 5)
        {
            tau->error_occurred = tmpl_True;
            tau->error_message = tmpl_strdup(
                "\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Get_Py_Perturb\n\n"
                "\rInput perturb is a list but does not have 5 entries.\n"
                "\rperturb must be a list of five real numbers.\n"
            );
            return;
        }

        iter = PyObject_GetIter(perturb);

        /*  Loop over the elements of the list, see if they can be converted  *
         *  to doubles, and store them in the tau->perturb variable.          */
        for (n = 0; n < 5; ++n)
        {
            next = PyIter_Next(iter);

            /*  If the element is an integer, convert to double and save it.  */
            if (PyLong_Check(next))
                tau->perturb[n] = PyLong_AsDouble(next);

            /*  Convert from Python float to C double with PyFloat_AsDouble.  */
            else if (PyFloat_Check(next))
                tau->perturb[n] = PyFloat_AsDouble(next);

            /*  Invalid data type for one of the entries. Return with error.  */
            else
            {
                tau->error_occurred = tmpl_True;
                tau->error_message = tmpl_strdup(
                    "\rError Encountered: rss_ringoccs\n"
                    "\r\trssringoccs_Get_Py_Perturb\n\n"
                    "\rInput perturb has entries that are not real numbers.\n"
                    "\rAll five entries for the perturb list must be numbers.\n"
                );
                return;
            }
        }
    }

    /*  The input was not a list. Return with error.                          */
    else
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Get_Py_Perturb\n\n"
            "\rInput perturb is not a list.\n"
        );
        return;
    }
}
