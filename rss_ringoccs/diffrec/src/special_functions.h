#ifndef RSS_RINGOCCS_SPECIAL_FUNCTIONS_H
#define RSS_RINGOCCS_SPECIAL_FUNCTIONS_H

/*  To avoid compiler warnings about deprecated numpy stuff.                  */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/*  complex data types, as well as _Complex_I, are defined here.              */
#include <complex.h>

/*  compute_norm_eq, max and min found here. math.h included here as well.    */
#include "__math_functions.h"
#include "__where.h"
#include "__get_array.h"
#include "__window_functions.h"

/*  Window functions and Fresnel transforms defined here.                     */
#include "__diffraction_functions.h"

/* All wrapper functions defined within these files.                          */
#include "_fraunhofer_diffraction_wrappers.h"
#include "__fresnel_diffraction.h"
#include "_fresnel_kernel_wrappers.h"
#include "_lambertw_wrappers.h"
#include "_physics_functions_wrappers.h"
#include "_resolution_inverse_function_wrappers.h"

/*  Various header files required for the C-Python API to work.               */
#include <Python.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

/*  Make sure these names are available.                                      */
#ifdef __get_one_real_from_one_real
#undef __get_one_real_from_one_real
#endif

#ifdef __get_one_real_from_two_real
#undef __get_one_real_from_two_real
#endif

#ifdef __get_one_real_from_three_real
#undef __get_one_real_from_three_real
#endif

#ifdef __get_complex_from_four_real
#undef __get_complex_from_four_real
#endif

#ifdef __get_complex_from_three_real
#undef __get_complex_from_three_real
#endif

/*  To avoid repeating the same code over and over again, define these macros *
 *  to be used for looping over functions.                                    */

#define __get_one_real_from_one_real(x, y, dim, f) ({\
    /*  Declare necessary variables.                                         */\
    long i;\
    \
    for (i=0; i<dim; ++i){\
        y[i] = (*f)(x[i]);\
    }\
})

#define __get_one_real_from_two_real(x1, x2, y, dim, f) ({\
    /*  Declare necessary variables.                                         */\
    long i;\
    \
    for (i=0; i<dim; ++i){\
        y[i] = (*f)(x1[i], x2);\
    }\
})

#define __get_one_real_from_three_real(x1, x2, x3, y, dim, f) ({\
    /*  Declare necessary variables.                                         */\
    long i;\
    \
    for (i=0; i<dim; ++i){\
        y[i] = (*f)(x1[i], x2, x3);\
    }\
})

#define __get_complex_from_three_real(x, a, F, y, dim, f) ({\
    /*  Declare necessary variables.                                         */\
    long i;\
    \
    for (i=0; i<dim; ++i){\
        y[i] = (*f)(x[i], a, F);\
    }\
})

#define __get_complex_from_four_real(x, a, b, F, y, dim, f) ({\
    /*  Declare necessary variables.                                         */\
    long i;\
    \
    for (i=0; i<dim; ++i){\
        y[i] = (*f)(x[i], a, b, F);\
    }\
})

#endif