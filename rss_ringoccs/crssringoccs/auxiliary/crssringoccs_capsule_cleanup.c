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

/*  Function prototype and typedefs for structs given here.                   */
#include "../crssringoccs.h"

/*  TMPL_FREE macro provided here.                                            */
#include <libtmpl/include/compat/tmpl_free.h>

/*  This function free's memory stored in certain Python objects, like numpy  *
 *  arrays. When the reference count for the object hits zero, this function  *
 *  is called and the buffer (allocated with malloc) is then free'd.          */
void crssringoccs_Capsule_Cleanup(PyObject * const capsule)
{
    void *memory = PyCapsule_GetPointer(capsule, NULL);
    TMPL_FREE(memory);
}
