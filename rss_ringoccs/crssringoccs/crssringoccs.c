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
#include "crssringoccs.h"

/*  Avoid warnings about deprecated Numpy API versions.                       */
#ifndef NPY_NO_DEPRECATED_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif

/*  Numpy header files.                                                       */
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    .m_name = "crssringoccs",
    .m_doc = "Module containing C Tools for rss_ringoccs.",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit_crssringoccs(void)
{
    PyObject *module = NULL;
    PyObject *all = NULL;
    int pymod_addobj;

    if (PyType_Ready(&crssringoccs_DiffractionCorrection) < 0)
        return NULL;

    if (PyType_Ready(&crssringoccs_ExtractCSVData) < 0)
        return NULL;

    if (PyType_Ready(&crssringoccs_GetUranusData) < 0)
        return NULL;

    if (PyType_Ready(&crssringoccs_GetMergedCSVData) < 0)
        return NULL;

    module = PyModule_Create(&moduledef);

    if (!module)
        return NULL;

    Py_INCREF(&crssringoccs_ExtractCSVData);
    pymod_addobj = PyModule_AddObject(
        module,
        "ExtractCSVData",
        (PyObject *)&crssringoccs_ExtractCSVData
    );

    if (pymod_addobj < 0)
    {
        Py_DECREF(&crssringoccs_ExtractCSVData);
        Py_DECREF(module);
        return NULL;
    }

    Py_INCREF(&crssringoccs_GetUranusData);
    pymod_addobj = PyModule_AddObject(
        module,
        "GetUranusData",
        (PyObject *)&crssringoccs_GetUranusData
    );

    if (pymod_addobj < 0)
    {
        Py_DECREF(&crssringoccs_ExtractCSVData);
        Py_DECREF(&crssringoccs_GetUranusData);
        Py_DECREF(module);
        return NULL;
    }

    Py_INCREF(&crssringoccs_DiffractionCorrection);
    pymod_addobj = PyModule_AddObject(
        module,
        "DiffractionCorrection",
        (PyObject *)&crssringoccs_DiffractionCorrection
    );

    if (pymod_addobj < 0)
    {
        Py_DECREF(&crssringoccs_ExtractCSVData);
        Py_DECREF(&crssringoccs_GetUranusData);
        Py_DECREF(&crssringoccs_DiffractionCorrection);
        Py_DECREF(module);
        return NULL;
    }

    Py_INCREF(&crssringoccs_GetMergedCSVData);
    pymod_addobj = PyModule_AddObject(
        module,
        "GetMergedCSVData",
        (PyObject *)&crssringoccs_GetMergedCSVData
    );

    if (pymod_addobj < 0)
    {
        Py_DECREF(&crssringoccs_ExtractCSVData);
        Py_DECREF(&crssringoccs_GetUranusData);
        Py_DECREF(&crssringoccs_DiffractionCorrection);
        Py_DECREF(&crssringoccs_GetMergedCSVData);
        Py_DECREF(module);
        return NULL;
    }

    all = Py_BuildValue(
        "[s, s, s, s]",
        "DiffractionCorrection",
        "ExtractCSVData",
        "GetUranusData",
        "GetMergedCSVData"
    );

    if (!all)
    {
        Py_DECREF(&crssringoccs_ExtractCSVData);
        Py_DECREF(&crssringoccs_GetUranusData);
        Py_DECREF(&crssringoccs_DiffractionCorrection);
        Py_DECREF(module);
        Py_CLEAR(module);
        return NULL;
    }

    pymod_addobj = PyModule_AddObject(module, "__all__", all);

    if (pymod_addobj < 0)
    {
        Py_DECREF(&crssringoccs_ExtractCSVData);
        Py_DECREF(&crssringoccs_GetUranusData);
        Py_DECREF(&crssringoccs_DiffractionCorrection);
        Py_DECREF(module);
        return NULL;
    }

    import_array();
    return module;
}
