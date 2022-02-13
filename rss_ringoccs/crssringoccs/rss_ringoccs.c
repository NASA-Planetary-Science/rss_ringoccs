/******************************************************************************
 *                                 LICENSE                                    *
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
#include <rss_ringoccs/include/crss_ringoccs.h>

/*  Numpy header files.                                                       */
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    .m_name = "custom",
    .m_doc = "Module containing C Tools for rss_ringoccs.",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit_crssringoccs(void)
{
    PyObject *m;
    int pymod_addobj;
    if (PyType_Ready(&DiffrecType) < 0)
        return NULL;

    if (PyType_Ready(&ExtractCSVDataType) < 0)
        return NULL;

    m = PyModule_Create(&moduledef);

    if (m == NULL)
        return NULL;

    Py_INCREF(&ExtractCSVDataType);
    pymod_addobj = PyModule_AddObject(m, "ExtractCSVData",
                                      (PyObject *) &ExtractCSVDataType);

    if (pymod_addobj < 0)
    {
        Py_DECREF(&ExtractCSVDataType);
        Py_DECREF(m);
        return NULL;
    }

    Py_INCREF(&DiffrecType);
    pymod_addobj = PyModule_AddObject(m, "DiffractionCorrection",
                                      (PyObject *) &DiffrecType);

    if (pymod_addobj < 0)
    {
        Py_DECREF(&DiffrecType);
        Py_DECREF(m);
        return NULL;
    }

    import_array();
    return m;
}

