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
#include <libtmpl/include/helper/tmpl_array_size.h>

/*  Avoid warnings about deprecated Numpy API versions.                       */
#ifndef NPY_NO_DEPRECATED_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif

/*  Tell numpy what the name of the C-Extension module is.                    */
#define PY_ARRAY_UNIQUE_SYMBOL crssringoccs

/*  Numpy header file containing the PyArray initialization routines.         */
#include <numpy/ndarrayobject.h>

/*  Struct containing a pointer to a C-Extension type (Python class) and the  *
 *  name of the type.                                                         */
struct crssringoccs_PythonType {
    PyTypeObject *type;
    const char *name;
};

/*  Array containing all of the C-Extension types. We'll loop over these      *
 *  entries and add them to the C-Extension module.                           */
struct crssringoccs_PythonType crssringoccs_TypeList[4] = {

    /*  DiffractionCorrection class, used for processing DLP data and         *
     *  creating reconstructed (diffraction corrected) profiles.              */
    {
        &crssringoccs_DiffractionCorrection,
        "DiffractionCorrection"
    },

    /*  CassiniCSVData class, reads GEO, CAL, and DLP .TAB files (CSV files)  *
     *  and loads them into a Python class that may be read by other classes  *
     *  such as DiffractionCorrection. These .TAB files must be in the PDS    *
     *  format (either the deprecated version or newer format) that is used   *
     *  for the Cassini data.                                                 */
    {
        &crssringoccs_CassiniCSVData,
        "CassiniCSVData"
    },

    /*  VoyagerCSVData class, similar to CassiniCSVData, reading in three     *
     *  files (GEO, CAL, DLP), but the DLP.TAB file must be in the format     *
     *  used for the Voyager data.                                            */
    {
        &crssringoccs_VoyagerCSVData,
        "VoyagerCSVData"
    },

    /*  MergedCSVData class, reads in a single Merged DLP file (DLPM.TAB) and *
     *  loads it into a Python class that may be used by other classes.       */
    {
        &crssringoccs_MergedCSVData,
        "MergedCSVData"
    }
};

/*  Definition of the C-Extension module, containg the name and docstring.    */
static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    .m_name = "crssringoccs",
    .m_doc = "Module containing C Tools for rss_ringoccs.",
    .m_size = -1,
};

/*  Helper function to add C-Extension types to the module.                   */
static int
crssringoccs_Add_Type(PyObject *module, struct crssringoccs_PythonType type)
{
    /*  Check if the type can be used by the Python interpreter.              */
    if (PyType_Ready(type.type) < 0)
        return 0;

    /*  The type is ready to use. Increment the reference counter to prevent  *
     *  it from being destroyed.                                              */
    Py_XINCREF(type.type);

    /*  Try to add the type to the module.                                    */
    if (PyModule_AddObject(module, type.name, (PyObject *)type.type) < 0)
    {
        /*  On error we need to decrement the reference counter to avoid      *
         *  possible memory leaks. Clear the type.                            */
        Py_CLEAR(type.type);
        return 0;
    }

    /*  Otherwise the type is ready and has been added to the module.         */
    return 1;
}

/*  Function for initializing the crssringoccs C-Extension module.            */
PyMODINIT_FUNC PyInit_crssringoccs(void)
{
    const size_t number_of_types = TMPL_ARRAY_SIZE(crssringoccs_TypeList);
    size_t m, n;
    int success;

    PyObject *all = NULL;
    PyObject *module = NULL;

    /*  Initialize numpy. This line is required before any use of the Numpy   *
     *  C API. segfaults are likely if this is omitted.                       */
    if (PyArray_ImportNumPyAPI() < 0)
        return NULL;

    /*  Create the C-Extension module.                                        */
    module = PyModule_Create(&moduledef);

    /*  PyModule_Create returns NULL on failure. Check for this.              */
    if (!module)
        return NULL;

    for (n = 0; n < number_of_types; ++n)
    {
        /*  Try to add the type to the module.                                */
        success = crssringoccs_Add_Type(module, crssringoccs_TypeList[n]);

        /*  Ensure that the type was successfully added.                      */
        if (!success)
        {
            /*  On error, clean up all previously added types.                */
            for (m = 0; m < n; ++m)
                Py_CLEAR(crssringoccs_TypeList[m].type);

            /*  Clear the module itself, it was not successfully created.     */
            Py_CLEAR(module);
            return NULL;
        }
    }

    /*  Expose all of the types to the user. These classes are accessible via *
     *  "from crssringoccs import ClassName". You may also use * to import    *
     *  everything using "from crssringoccs import *"".                       */
    all = Py_BuildValue(
        "[s, s, s, s]",
        "DiffractionCorrection",
        "CassiniCSVData",
        "VoyagerCSVData",
        "MergedCSVData"
    );

    /*  Py_BuildValue returns NULL on failure. Check for this.                */
    if (!all)
    {
        /*  On error, clean up everything that was initialized.               */
        for (n = 0; n < number_of_types; ++n)
            Py_CLEAR(crssringoccs_TypeList[n].type);

        /*  Lastly, clean up the module and return with error.                */
        Py_CLEAR(module);
        return NULL;
    }

    /*  Otherwise, add this "all" object as the __all__ attribute for the     *
     *  module. This allows the C types to be easily accessed.                */
    if (PyModule_AddObject(module, "__all__", all) < 0)
    {
        /*  On error, clear all of the types that were initialized.           */
        for (n = 0; n < number_of_types; ++n)
            Py_CLEAR(crssringoccs_TypeList[n].type);

        /*  Clean up the "all" object and the module as well.                 */
        Py_CLEAR(all);
        Py_CLEAR(module);
        return NULL;
    }

    return module;
}
/*  End of PyInit_crssringoccs.                                               */
