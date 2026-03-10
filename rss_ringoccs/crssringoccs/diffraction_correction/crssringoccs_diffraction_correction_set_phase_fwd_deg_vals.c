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
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       March 10, 2026                                                *
 ******************************************************************************/

/*  Function prototype and typedefs for structs given here.                   */
#include "../crssringoccs.h"

int
crssringoccs_DiffractionCorrection_Set_Phase_Fwd_Deg_Vals(PyObject *op,
                                                          PyObject *value,
                                                          void *closure)
{
    /*  Variable for the diffraction correction class instance.               */
    crssringoccs_PyDiffrecObj *self;

    /*  A NULL pointer for the first argument should be treated as an error.  */
    if (!op)
    {
        PyErr_Format(
            PyExc_RuntimeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_DiffractionCorrection_Set_Phase_Fwd_Deg_Vals\n\n"
            "\rInput 'op' is NULL.\n\n"
        );

        return -1;
    }

    /*  Get a pointer to the actual DiffractionCorrection instance.           */
    self = (crssringoccs_PyDiffrecObj *)op;

    /*  Set phase_fwd_deg_vals to the input value. There are no checks for    *
     *  what type of object this is, the user is free to set it to anything.  */
    Py_XSETREF(self->phase_fwd_deg_vals, Py_XNewRef(value));
    return 0;
}
/*  End of crssringoccs_DiffractionCorrection_Set_Phase_Fwd_Deg_Vals.         */
