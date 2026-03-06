
/*  Include guard to prevent including this file twice.                       */
#ifndef CRSSRINGOCCS_NUMPY_API_H
#define CRSSRINGOCCS_NUMPY_API_H

/*  Avoid warnings about deprecated Numpy API versions. See:                  *
 *      https://numpy.org/doc/stable/reference/c-api/deprecations.html        *
 *  for more details.                                                         */
#ifndef NPY_NO_DEPRECATED_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif

/*  Define both of the macros so that the helper functions do not need to     *
 *  repeatedly call NumPy's import_array macro. See:                          *
 *      https://numpy.org/devdocs/reference/c-api/array.html                  *
 *  for details. The name of the PY_ARRAY_UNIQUE_SYMBOL macro is the name of  *
 *  the compiled Python C-Extension module, which is crssringoccs, and not    *
 *  the standalone C library (librssringoccs).                                *
 *                                                                            *
 *  Note, this file should not be included into the main module file,         *
 *  crssringoccs.c. That file should not have NO_IMPORT_ARRAY defined at all. */
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL crssringoccs

/*  PyArrayObject typedef and the Numpy-C Array API provided here.            */
#include <numpy/ndarraytypes.h>
#include <numpy/ndarrayobject.h>

#endif
/*  End of include guard.                                                     */
