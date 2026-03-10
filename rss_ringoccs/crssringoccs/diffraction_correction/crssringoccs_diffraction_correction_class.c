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

PyTypeObject crssringoccs_DiffractionCorrection = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "DiffractionCorrection",
    .tp_doc = "Diffraction Correction class.",
    .tp_basicsize = sizeof(crssringoccs_PyDiffrecObj),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = crssringoccs_DiffractionCorrection_New,
    .tp_init = (initproc)crssringoccs_DiffractionCorrection_Init,
    .tp_dealloc = (destructor)crssringoccs_DiffractionCorrection_Destroy,
    .tp_getset = crssringoccs_DiffractionCorrection_GetSetters,
    .tp_members = crssringoccs_DiffractionCorrection_Members,
    .tp_methods = crssringoccs_DiffractionCorrection_Methods
};
