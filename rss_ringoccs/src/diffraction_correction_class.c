/******************************************************************************
 *                                 LICENSE                                    *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify it   *
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
 *                             Diffrec Module                                 *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      This file contains the C-Python API code needed to use the C code     *
 *      contained in the diffrec/ directory within Python. It contains        *
 *      wrappers for several functions and the DiffractionCorrection class,   *
 *      which is the primary class for performing diffraction correction.     *
 *      Upon successful compilation there should be a diffrec.*info*.so file. *
 *      The Python interpreter can directly import diffrec if this file is in *
 *      its path via:                                                         *
 *          import diffrec                                                    *
 *      More likely, the user will be at the top-level of rss_ringoccs, where *
 *      one can simply run the following:                                     *
 *          from rss_ringoccs import diffrec                                  *
 *      The following are defined in this file:                               *
 *          DiffractionCorrection:                                            *
 *              A Python class, written in C using the C-Python API. The      *
 *              equivalent of the __init__ function at the python level is    *
 *              defined via Diffrec_init. This is the main class for          *
 *              processing diffracted data. See the file                      *
 *              _diffraction_correction.c and the associated header file      *
 *              _diffraction_correction.h for more information on the         *
 *              algorithms implemented.                                       *
 *          gap_diffraction:                                                  *
 *              A diffraction modeling tool which returns the diffraction     *
 *              profile of a gap, which has 0 complex amplitude everywhere    *
 *              except the gap where it evaluates to 1.                       *
 *          ringlet_diffraction:                                              *
 *              Similar to gap diffraction with the values reversed.          *
 *          right_straightedge:                                               *
 *              diffraction modeling of a "straight-edge". This has value     *
 *              1 for everything to the right of the straight edge, and 0 to  *
 *              the left.                                                     *
 *          left_straightedge:                                                *
 *              Like right_straightedge with the values reverse.              *
 *          fresnel_psi:                                                      *
 *              The Fresnel kernel function, as defined in MTR86.             *
 *          fresnel_scale:                                                    *
 *              The Fresnel scale, also defined in MTR86.                     *
 *          fresnel_dpsi_dphi:                                                *
 *              The partial derivative of the Fresnel kernel with respect to  *
 *              the phi variable.                                             *
 *          fresnel_dphi_dphi_ellipse:                                        *
 *              The partial derivative of the Fresnel kernel with respect to  *
 *              the phi variable where the kernel has been perturbed by       *
 *              the elliptical geometry of the rings.                         *
 *          fresnel_d2phi_dphi2:                                              *
 *              The second partial derivative of the Fresnel kernel with      *
 *              respect to the phi variable.                                  *
 *          single_slit_diffraction:                                          *
 *              Fraunhofer diffraction modeling of a single slit.             *
 *          double_slit_diffraction:                                          *
 *              Fraunhofer diffraction modeling of a double slit.             *
 *          square_wave_diffraction:                                          *
 *              Fresnel diffraction modeling of a square wave.                *
 *          square_well_phase:                                                *
 *              The modeled phase of the diffraction profile from a ringlet.  *
 ******************************************************************************
 *                              DEFINED FUNCTIONS                             *
 ******************************************************************************
 *  The following functions are defined and used throughout, but are not      *
 *  available to the Python interpreter.                                      *
 *  GetNormeqFromString:                                                      *
 *      This acts like a Python dictionary where we input a string and        *
 *      retrieve a predefined value. It gives us the normalized equivalent    *
 *      width corresponding to a legal string.                                *
 *      NOTE:                                                                 *
 *          If you want to add your own window function, and your own         *
 *          corresponding normalized equivalent width, you will also need to  *
 *          add the corresponding string (the name of the window) to several  *
 *          of the error checking routines throughout the code.               *
 *  GetRangeFromString:                                                       *
 *      This also acts like a Python dictionary, giving us two values, the    *
 *      start and end points, corresponding to a given legal string.          *
 *      NOTE:                                                                 *
 *          Currently, all legal strings correspond to features of Saturn.    *
 *          If you wish to add your own string/feature, you will need to      *
 *          append it to this function and the error checking routines.       *
 *  capsule_cleanup:                                                          *
 *      A very important function, even though it is only two lines. This     *
 *      function ensures that memory we have allocated for variables in C     *
 *      that are then passed to the Python interpreter will still be freed    *
 *      when Python no longer needs it. Without this we will have serious     *
 *      memory leaks. It is primiarly used for safely returning numpy arrays. *
 ******************************************************************************
 *                            A FRIENDLY WARNING                              *
 ******************************************************************************
 *  All of the functions are declared "static" and cannot be used outside of  *
 *  this file. The only purpose of this file is generating the diffrec module.*
 *                                                                            *
 *  This code uses complex numbers throughout, and is compatible with the C99 *
 *  standard. To use this code, make sure your compiler supports C99 or more  *
 *  recent standards of the C Programming Language.                           *
 *                                                                            *
 *  In addition, for full use of all the routines one will need the FFTW      *
 *  library. Installation instruction can be found in the rss_ringoccs PDF.   *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  For successful compiling you will need the following NON-STANDARD         *
 *  header files / libraries:                                                 *
 *  libfftw3:                                                                 *
 *      The shared library file for the FFTW package. This is needed to       *
 *      use the FFT options of the DiffractionCorrection class.               *
 *  Python.h:                                                                 *
 *      The main header file for Python. If this code is compiled via the     *
 *      setup.py script provided, then Python.h will automatically be         *
 *      included correctly.                                                   *
 *  numpy/ndarraytypes.h:                                                     *
 *  numpy/ufuncobject.h:                                                      *
 *      Header files that allow us to use the Numpy API functions for         *
 *      accessing and creating numpy arrays at the C level. To get the        *
 *      location of these files you will need to use numpy's get_include()    *
 *      function. The setup.py script does this for you.                      *
 *  librssringoccsdiffractioncorrection.so:                                   *
 *      A non-standard library written and maintained by the rss_ringoccs     *
 *      project, it provide all of the diffraction correction algorithms      *
 *      written at the C level. The setup.py file will automatically link     *
 *      this program to the correct headers. To use this module Python will   *
 *      need to be able to find this .so. The setup scripts will install it   *
 *      in /usr/local/lib which is a standard location on Linux and MacOS.    *
 *  librssringoccsspecialfunctions.so:                                        *
 *      Another non-standard library maintained by rss_ringoccs, it contains  *
 *      a plethora of mathematical functions that are used in diffraction     *
 *      correction such as bessel functions, Legendre polynomials, and many   *
 *      more. The same criterion for librssringoccsdiffractioncorrection must *
 *      be satisfied in order to correctly link and load this library.        *
 *      This file is required implicitly by the                               *
 *      librssringoccsdiffractioncorrection.so library. We do not directly    *
 *      make calls to its functions.                                          *
 ******************************************************************************
 *                            A NOTE ON COMMENTS                              *
 ******************************************************************************
 *  It is anticipated that many users of this code will have experience in    *
 *  either Python or IDL, but not C. Many comments are left to explain as     *
 *  much as possible. Vagueness or unclear code should be reported to:        *
 *  https://github.com/NASA-Planetary-Science/rss_ringoccs/issues             *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       June 22, 2019                                                 *
 ******************************************************************************
 *                                History                                     *
 ******************************************************************************
 *  2020/09/01 (Ryan Maguire):                                                *
 *      Added DiffractionCorrection class in C. Updated comments.             *
 ******************************************************************************/

/*  To avoid compiler warnings about deprecated numpy stuff.                  */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/*  The standard library header stdlib contains malloc, calloc, and realloc,  *
 *  free, as well as strtol (string-to-long) which converts a string like     *
 *  "314" to the integer 314.                                                 */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

/*  To help manipulate strings, both native C strings const char* and char*,  *
 *  as well as Python strings passed by the user, include string.h. This is a *
 *  standard library header containing strcmp (string-compare), strncmp       *
 *  (compares first n characters of string), strcpy (string-copy), strncpy,   *
 *  and strlen (string-length), all of which are used at some point.          */
#include <string.h>

/*  Useful library for manipulate strings, it contains the tolower function   *
 *  which takes a string like "ABC123abc" and return "abc123abc". Since all   *
 *  of the strings passed to the diffraction correction routines need to be   *
 *  lower-case, we will make use of this frequently.                          */
#include <ctype.h>

/*  The following are NON-STANDARD header files that MUST BE IN YOUR PATH.    *
 *  If you installed python using anaconda then Python.h should automatically *
 *  be included in your path. Also, if you are using the setup.py script      *
 *  provided then inclusion of these files should be done for you.            */
#include <Python.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

/*  The following header files are NON-STANDARD, and are a part of the        *
 *  rss_ringoccs package. The setup scripts will add the correct CFLAGS so    *
 *  compiler should see these without user interaction.                       */

/*  This contains the PyDiffrecObj, using the C-Python API, which allows us   *
 *  to define the DiffractionCorrection class in C.                           */
#include "rss_ringoccs_py_api.h"
#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>

/*  To avoid repeating the same code over and over again, we'll define a few  *
 *  preprocessor functions (macros). First make sure the names are available. */

/*  This macro takes a the word Var and converts to a string at compile time. */
#ifdef VarToString
#undef VarToString
#endif

#ifdef _array_from_three
#undef _array_from_three
#endif

#ifdef _array_from_four
#undef _array_from_four
#endif

#define VarToString(Var) (#Var)

#define _array_from_three(x1, x2, x3, y, dim, intype, outtype, f) ({           \
    unsgiend long i;                                                           \
    outtype *out = (outtype *)malloc(sizeof(outtype)*dim);                     \
    for (i=0; i<dim; ++i)                                                      \
        out[i] = (*f)(((intype *)x1)[i], x2, x3);                              \
    y = out;                                                                   \
})

#define _array_from_four(x1, x2, x3, x4, y, dim, intype, outtype, f) ({        \
    unsigned long i;                                                           \
    outtype *out = (outtype *)malloc(sizeof(outtype)*dim);                     \
    for (i=0; i<dim; ++i)                                                      \
        out[i] = (*f)(((intype *)x1)[i], x2, x3, x4);                          \
    y = out;                                                                   \
})

/*  For when a value/array so be non-negative.                                */
#define RSS_RINGOCCSNegValueError(FuncName, VarName)                           \
PyErr_Format(                                                                  \
    PyExc_ValueError,                                                          \
    "\n\rError Encountered: rss_ringoccs\n"                                    \
    "\r\tdiffrec.%s\n\n"                                                       \
    "\r%s must be positive.",                                                  \
    VarToString(FuncName), VarToString(VarName)                                \
);

/*  Raise this when a value should be between -2pi and 2pi.                   */
#define RSS_RINGOCCSTwoPiValueError(FuncName, VarName)                         \
PyErr_Format(                                                                  \
    PyExc_ValueError,                                                          \
    "\n\rError Encountered: rss_ringoccs\n"                                    \
    "\r\tdiffrec.%s\n\n"                                                       \
    "\r%s must be between minus two pi and two pi (radians).",                 \
    VarToString(FuncName), VarToString(VarName)                                \
);

/*  Macro for raising the appropriate python error if the DLP instance is     *
 *  missing an attribute. This is equivalent to the following in python       *
 *      if not hasattr(tauin, attr_name):                                     *
 *          raise AttributeError(                                             *
 *              """                                                           *
 *              Error message                                                 *
 *              """                                                           *
 *          )                                                                 *
 *      else:                                                                 *
 *          pass                                                              *
 *  It then checks that the variable is a numpy array using numpy's API. This *
 *  is equivalent to the following:                                           *
 *      if not isinstance(varname, numpy.ndarray):                            *
 *          raise TypeError(                                                  *
 *              """                                                           *
 *              Error message                                                 *
 *              """                                                           *
 *          )                                                                 *
 *      else:                                                                 *
 *          pass                                                              *
 *  Next we try to convert the numpy array to an array of double, which is    *
 *  equivalent to using the astype method of the ndarray numpy object:        *
 *      arr = arr.astype(float)                                               *
 *  Finally, we check that the array is one dimensional and that it has the   *
 *  same number of elements as the input rho_km_vals array. If this passes,   *
 *  we pointer the pointer ptr to the data of the array.                      */
#define RSS_RINGOCCSCheckArray(ptr, tauin, attr_name, tmp, arr, length)        \
if (!PyObject_HasAttrString(tauin, VarToString(attr_name)))                    \
{                                                                              \
    PyErr_Format(                                                              \
        PyExc_AttributeError,                                                  \
        "\n\rError Encountered: rss_ringoccs\n"                                \
        "\r\tdiffrec.DiffractionCorrection\n\n"                                \
        "\rInput DLP Instance is missing the following attribute:\n"           \
        "\r\t%s\n\n", VarToString(attr_name)                                   \
    );                                                                         \
    return -1;                                                                 \
}                                                                              \
else                                                                           \
    tmp = PyObject_GetAttrString(tauin, VarToString(attr_name));               \
                                                                               \
if (!PyArray_Check(tmp))                                                       \
{                                                                              \
    PyErr_Format(                                                              \
        PyExc_TypeError,                                                       \
        "\n\rError Encountered: rss_ringoccs\n"                                \
        "\r\tdiffrec.DiffractionCorrection\n\n"                                \
        "\r%s must be a numpy array.\n",                                       \
        VarToString(varname)                                                   \
    );                                                                         \
    return -1;                                                                 \
}                                                                              \
else                                                                           \
    arr = PyArray_FromObject(tmp, NPY_DOUBLE, 1, 1);                           \
                                                                               \
/*  If PyArray_FromObject failed arr should be NULL. If so, raise error.     */\
if (!arr)                                                                      \
{                                                                              \
    PyErr_Format(                                                              \
        PyExc_TypeError,                                                       \
        "\n\rError Encountered: rss_ringoccs\n"                                \
        "\r\tdiffrec.DiffractionCorrection\n\n"                                \
        "\rCould not convert %s to double array. Input is most\n",             \
        "\rlikely complex numbers or contains a string.\n\n",                  \
        VarToString(varname)                                                   \
    );                                                                         \
    return -1;                                                                 \
}                                                                              \
                                                                               \
/*  Currently we only allow for one dimensional inputs.                      */\
else if (PyArray_NDIM((PyArrayObject *)arr) != 1)                              \
{                                                                              \
    PyErr_Format(                                                              \
        PyExc_IndexError,                                                      \
        "\n\rError Encountered: rss_ringoccs\n"                                \
        "\r\tdiffrec.DiffractionCorrection\n\n"                                \
        "\r%s must be a one-dimensional numpy array.\n",                       \
        VarToString(varname)                                                   \
    );                                                                         \
    return -1;                                                                 \
}                                                                              \
                                                                               \
/*  arr should have the same number of elements as rho_km_vals, which we are */\
/*  passing to this preprocessor functions as "length".                      */\
else if (PyArray_DIMS((PyArrayObject *)arr)[0] != (long)length)                \
{                                                                              \
    PyErr_Format(                                                              \
        PyExc_IndexError,                                                      \
        "\n\rError Encountered: rss_ringoccs\n"                                \
        "\r\tdiffrec.DiffractionCorrection\n\n"                                \
        "\r%s and rho_km_vals have a different number of elements.\n",         \
        VarToString(varname)                                                   \
    );                                                                         \
    return -1;                                                                 \
}                                                                              \
                                                                               \
/*  If every passed, set ptr to point to the data inside the array arr.      */\
else                                                                           \
    ptr = (double *)PyArray_DATA((PyArrayObject *)arr);

/*  Use this for setting the data into the DiffractionCorrection class.      */\
#define RSS_RINGOCCSSetArray(self, tau, VarName, tmp, typenum, start, length)  \
PyObject *PyArray_##VarName;                                                   \
PyObject *PyCapsule_##VarName;                                                 \
double *VarName = malloc(sizeof(double) * length);                             \
long pylength = (long)length;                                                  \
                                                                               \
for (i=0; i<length; ++i)                                                       \
    VarName[i] = tau.VarName[i+start];                                         \
                                                                               \
PyArray_##VarName = PyArray_SimpleNewFromData(1, &pylength, typenum,           \
                                             (void *)VarName);                 \
PyCapsule_##VarName = PyCapsule_New((void *)VarName, NULL, capsule_cleanup);   \
                                                                               \
/*  This frees the variable at the Python level once it's destroyed.         */\
PyArray_SetBaseObject((PyArrayObject *)PyArray_##VarName, PyCapsule_##VarName);\
                                                                               \
tmp = self->VarName;                                                           \
Py_INCREF(PyArray_##VarName);                                                  \
self->VarName = PyArray_##VarName;                                             \
Py_XDECREF(tmp);

/*  Macro for checking if all of the elements of an array are positive.       */
#define RSS_RINGOCCSCheckArrayPositive(varname, arrsize)                       \
if (rssringoccs_Min_Double(varname, arrsize) <= 0.0)                           \
{                                                                              \
    PyErr_Format(                                                              \
        PyExc_IndexError,                                                      \
        "\n\rError Encountered: rss_ringoccs\n"                                \
        "\r\tdiffrec.DiffractionCorrection\n\n"                                \
        "\r%s contains non-positive values.\n",                                \
        VarToString(varname)                                                   \
    );                                                                         \
    return -1;                                                                 \
}

/*  Macro for checking if the values of a numpy array fall between -2pi,2pi.  */
#define RSS_RINGOCCSCheckArrayTwoPi(varname, arrsize)                          \
if ((rssringoccs_Min_Double(varname, arrsize) < -rssringoccs_Two_Pi) ||        \
    (rssringoccs_Max_Double(varname, arrsize) > rssringoccs_Two_Pi))           \
{                                                                              \
    PyErr_Format(                                                              \
        PyExc_IndexError,                                                      \
        "\n\rError Encountered: rss_ringoccs\n"                                \
        "\r\tdiffrec.DiffractionCorrection\n\n"                                \
        "\r%s contains values outside of [-2pi, 2pi].\n",                      \
        VarToString(varname)                                                   \
    );                                                                         \
    return -1;                                                                 \
}

#define RSS_RINGOCCSRangeError(range_input)                                    \
    PyErr_Format(                                                              \
        PyExc_ValueError,                                                      \
        "\n\rError Encountered: rss_ringoccs\n"                                \
        "\r\tdiffrec.DiffractionCorrection\n\n"                                \
        "\rIllegal string for rng. Allowed strings are:\n"                     \
        "\r\tall               [1.0, 400000.0]\n"                              \
        "\r\tbesselbarnard     [120210.0, 120330.0]\n"                         \
        "\r\tbessel-barnard    [120210.0, 120330.0]\n"                         \
        "\r\tcringripples      [77690.0, 77760.0]\n"                           \
        "\r\tencke             [132900.0, 134200.0]\n"                         \
        "\r\tenckegap          [132900.0, 134200.0]\n"                         \
        "\r\therschel          [118100.0, 118380.0]\n"                         \
        "\r\therschelgap       [118100.0, 118380.0]\n"                         \
        "\r\thuygens           [117650.0, 117950.0]\n"                         \
        "\r\thuygensringlet    [117650.0, 117950.0]\n"                         \
        "\r\tjanusepimetheus   [96200.0, 96800.0]\n"                           \
        "\r\tjeffreys          [118900.0, 119000.0]\n"                         \
        "\r\tjeffreysgap       [118900.0, 119000.0]\n"                         \
        "\r\tkuiper            [119300.0, 119500.0]\n"                         \
        "\r\tkuipergap         [119300.0, 119500.0]\n"                         \
        "\r\tmaxwell           [87410.0, 87610.0]\n"                           \
        "\r\tmaxwellringlet    [87410.0, 87610.0]\n"                           \
        "\r\trussell           [118550.0, 118660.0]\n"                         \
        "\r\trussellgap        [118550.0, 118660.0]\n"                         \
        "\r\ttitan             [77870.0, 77930.0]\n"                           \
        "\r\ttitanringlet      [77870.0, 77930.0\n\n"                          \
        "\rYour string:         %s",                                           \
        range_input                                                            \
    );

/*  This function frees the memory allocated to a pointer by malloc when the  *
 *  corresponding variable is destroyed at the Python level. Without this you *
 *  will have serious memory leaks, so do not remove!                         */
static void capsule_cleanup(PyObject *capsule)
{
    void *memory = PyCapsule_GetPointer(capsule, NULL);
    free(memory);
}

/*  Deallocating function for the DiffractionCorrection class.                */
static void Diffrec_dealloc(PyDiffrecObj *self)
{
    Py_XDECREF(self->B_rad_vals);
    Py_XDECREF(self->D_km_vals);
    Py_XDECREF(self->F_km_vals);
    Py_XDECREF(self->f_sky_hz_vals);
    Py_XDECREF(self->outfiles);
    Py_XDECREF(self->p_norm_fwd_vals);
    Py_XDECREF(self->p_norm_vals);
    Py_XDECREF(self->phase_fwd_vals);
    Py_XDECREF(self->phase_rad_vals);
    Py_XDECREF(self->phase_vals);
    Py_XDECREF(self->phi_rad_vals);
    Py_XDECREF(self->phi_rl_rad_vals);
    Py_XDECREF(self->power_vals);
    Py_XDECREF(self->raw_tau_threshold_vals);
    Py_XDECREF(self->rev_info);
    Py_XDECREF(self->rho_corr_pole_km_vals);
    Py_XDECREF(self->rho_corr_timing_km_vals);
    Py_XDECREF(self->rho_dot_kms_vals);
    Py_XDECREF(self->rho_km_vals);
    Py_XDECREF(self->t_oet_spm_vals);
    Py_XDECREF(self->t_ret_spm_vals);
    Py_XDECREF(self->t_set_spm_vals);
    Py_XDECREF(self->tau_threshold_vals);
    Py_XDECREF(self->tau_vals);
    Py_XDECREF(self->w_km_vals);
    Py_XDECREF(self->dathist);
    Py_XDECREF(self->history);
    Py_XDECREF(self->perturb);
    Py_XDECREF(self->rngreq);
    Py_TYPE(self)->tp_free((PyObject *) self);
}

/*  The init function for the dirrection correction class. This is the        *
 *  equivalent of the __init__ function defined in a normal python class.     */
static int Diffrec_init(PyDiffrecObj *self, PyObject *args, PyObject *kwds)
{
    /*  The list of the keywords accepted by the DiffractionCorrection class. *
     *  dlp and res are REQUIRED inputs, the rest are optional. If the user   *
     *  does not provide these optional keywords, we must set them ourselves. */
    static char *kwlist[] = {
        "dlp",
        "res",
        "rng",
        "wtype",
        "fwd",
        "use_norm",
        "use_fft",
        "verbose",
        "bfac",
        "sigma",
        "psitype",
        "write_file",
        "res_factor",
        "ecc",
        "peri",
        "perturb",
        "interp",
        NULL
    };

    /*  Variable for indexing.                                                */
    unsigned long i;

    /*  Python objects needed throughout the computation.                     */
    PyObject *DLPInst;
    PyObject *iter;
    PyObject *next;
    PyObject *tmp;
    PyObject *arr;
    PyObject *history;
    PyObject *rngreq;
    PyObject *perturb;

    /*  Set the default keyword options.                                      */

    /*  The kbmd20 is a new window, a modifed Kaiser-Bessel with alpha set to *
     *  two pi. The modification ensures the window goes to zero at its edges *
     *  while evaluating to one at the center, unlike the actual              *
     *  Kaiser-Bessel which is discontinuous at the edge of the window. The   *
     *  two pi factor is mostly guess work since it accurately reproduces the *
     *  PDS results. The real window used for that data is not known too me.  *
     *  The actual code for the window functions is in special_functions/     */
    self->wtype = "kbmd20";

    /*  Fresnel 4 is a new option, not mentioned in any of the papers but     *
     *  documented in our accompanying PDF. It uses Legendre polynomials to   *
     *  approximate the Fresnel kernel. It essentially takes Fresnels         *
     *  quadratic method to the next step, a quartic, hence the name. It is   *
     *  extremely fast (all of Rev007 takes less than a second) and very      *
     *  accurate for all but the most extreme occultations (like Rev133).     */
    self->psitype = "fresnel4";

    /*  Default range is "all", denoting [1.0, 400000.0]. We'll set later.    */
    rngreq = NULL;

    /*  By default, forward computations are not run, FFTs are not used, and  *
     *  the run is silent (verbose is off).                                   */
    self->fwd = rssringoccs_False;
    self->use_fft = rssringoccs_False;
    self->verbose = rssringoccs_False;

    /*  Using the bfac guarantees accurate window sizes in the case of a poor *
     *  Allen deviation. Window normalization is also recommended since the   *
     *  integral is scaled by the width of the window, and hence for small    *
     *  window sizes the result might return close to zero.                   */
    self->bfac = rssringoccs_True;
    self->use_norm = rssringoccs_True;

    /*  The default sigma value is the one for Cassini.                       */
    self->sigma = 2.0e-13;

    /*  If res_factor was not set, set to 0.75. This value was specified by   *
     *  Essam Marouf as necessary to ensure the reconstruction matches the    *
     *  PDS results. No justification is known to me.                         */
    self->res_factor = 0.75;

    /*  The write files option allows the user to export the end result to a  *
     *  CSV. Default is set to false.                                         */
    self->write_file = rssringoccs_False;

    /*  The default geometry assumes the rings are circular, so we set both   *
     *  the eccentricity and the periapse to zero.                            */
    self->ecc = 0.0;
    self->peri = 0.0;

    /*  By default, polynomial interpolation is off.                          */
    self->interp = 0;

    /*  Default polynomial perturbation is off.                               */
    perturb = NULL;

    /*  The perturb variable is a list and should contain exactly 5 elements. */
    const int PERTURB_LEN = 5;

    /*  If the range input is a list, it should have exactly two elements.    */
    const int RANGE_LEN = 2;

    /*  Useful numbers for checking the input psitype and wtype.              */
    const int FRESNEL_STR_LENGTH = 7;
    const int BASETEN = 10;

    /*  Extract the inputs and keywords supplied by the user. If the data     *
     *  cannot be extracted, raise a type error and return to caller. A short *
     *  explaination of PyArg_ParseTupleAndKeywords. The inputs args and kwds *
     *  are somewhat straight-forward, they're the arguments and keywords     *
     *  passed by the string. The cryptic string is not straight-forward. The *
     *  | symbol means everything after need not be positional, and we can    *
     *  specify arguments and keywords by name when calling                   *
     *  DiffractionCorrection, for example                                    *
     *  DiffractionCorrect(..., wtype="blah"). O indicates a Python object,   *
     *  and d is a Python float. This is the DLP and res variables. The $     *
     *  symbold means everything after is optional. s is a string, p is a     *
     *  Boolean (p for "predicate"). b is an integer, and the colon : denotes *
     *  that the input list has ended.                                        */
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|Od$OspppppdspdddOb:", kwlist,
                                     &DLPInst,          &self->input_res,
                                     &rngreq,           &self->wtype,
                                     &self->fwd,        &self->use_norm,
                                     &self->use_fft,    &self->verbose,
                                     &self->bfac,       &self->sigma,
                                     &self->psitype,    &self->write_file,
                                     &self->res_factor, &self->ecc,
                                     &self->peri,       &perturb,
                                     &self->interp))
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rCould not parse input variables.\n\n"
            "\rInputs:\n"
            "\r\tDLPInst:  \tAn instance of the DLP Class.\n"
            "\r\tres:      \tRequested resolution in km (float).\n\n"
            "\rKeywords:\n"
            "\rrng         \tThe requested range (str or list).\n"
            "\rwtype       \tThe requested window type (str).\n"
            "\rfwd         \tForward computation (bool).\n"
            "\ruse_norm    \tWindow normalization (bool).\n"
            "\ruse_fft     \tComputation using FFTs (bool).\n"
            "\rverbose     \tPrint status updates (bool).\n"
            "\rbfac        \tUse b-factor in window width (bool).\n"
            "\rsigma       \tThe Allen deviation (float).\n"
            "\rpsitype     \tRequested Frensel kernel approximation (str).\n"
            "\rwrite_file  \tWrite output to file (bool).\n"
            "\rres_factor  \tScaling factor for resolution (float).\n"
            "\recc         \tEccentricity of rings (bool).\n"
            "\rperi        \tPeriapse of rings (bool).\n"
            "\rperturb     \tRequested perturbation to Fresnel kernel (list).\n"
            "\rinterp      \tPolynomial interpolation to kernel (int).\n"
        );
        return -1;
    }

    /*  If verbose was set, print a status update.                            */
    if (self->verbose)
    {
        puts("Processing Diffraction Correction");
        puts("\tRunning Error Check on Input Arguments...");
    }

    /*  Create an instance of the DLPObj (C struct).                          */
    rssringoccs_TAUObj tau;

    /*  Check that interp is a legal value.                                   */
    if ((0 <= self->interp) && (self->interp <=4))
        tau.interp = self->interp;
    else
    {
        PyErr_Format(
            PyExc_ValueError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rinterp should be an integer between 0 and 4.\n"
            "\rYour input: %d\n", self->interp
        );
        return -1;
    }

    /*  Check that the input perturb is a list with 5 elements.               */
    if (!(perturb))
    {
        /*  No perturb was provided by the user, so set everything to None.   */
        tmp = self->perturb;
        perturb = Py_None;

        /*  tau.perturb still needs to be set, so loop over and set to zero.  */
        for (i=0; i<PERTURB_LEN; ++i)
            tau.perturb[i] = 0.0;

        /*  Safely store perturb in self. This method is recommended in the   *
         *  C-Python API documentation.                                       */
        Py_INCREF(perturb);
        self->perturb = perturb;
        Py_XDECREF(tmp);
    }

    /*  If the user supplied a perturb list, parse it and extract values.     */
    else if (PyList_Check(perturb))
    {
        /*  If the list is not the correct size, raise an error.              */
        if (!(PyList_Size(perturb) == PERTURB_LEN))
        {
            PyErr_Format(
                PyExc_IndexError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tdiffrec.DiffractionCorrection\n\n"
                "\rperturb should be a python list of five floats or ints.\n"
                "\rYour list has %d elements\n", PyList_Size(self->perturb)
            );
            return -1;
        }

        iter = PyObject_GetIter(perturb);
        if (!iter)
        {
            /*  The list is empty, return with an error.                      */
            PyErr_Format(
                PyExc_IndexError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tdiffrec.DiffractionCorrection\n\n"
                "\rperturb should be a python list of five floats or ints.\n"
            );
            return -1;
        }

        /*  Loop over the elements of the list, see if they can be converted  *
         *  to doubles, and store them in the tau.perturb variable.           */
        for (i=0; i<PERTURB_LEN; ++i)
        {
            next = PyIter_Next(iter);
            if (!next)
            {
                /*  The list is not big enough, return with error.            */
                PyErr_Format(
                    PyExc_IndexError,
                    "\n\rError Encountered: rss_ringoccs\n"
                    "\r\tdiffrec.DiffractionCorrection\n\n"
                    "\rperturb should be a python list of five floats or\n"
                    "\rints. Your list has size: %d", i
                );
                return -1;
            }

            /*  If the element is an integer, convert to double and save it.  */
            if (PyLong_Check(next))
                tau.perturb[i] = (double)PyLong_AsLong(next);

            /*  Convert from Python float to C double with PyFloat_AsDouble.  */
            else if (PyFloat_Check(next))
                tau.perturb[i] = PyFloat_AsDouble(next);

            /*  Invalid data type for one of the entries. Return with error.  */
            else
            {
                PyErr_Format(
                    PyExc_TypeError,
                    "\n\rError Encountered: rss_ringoccs\n"
                    "\r\tdiffrec.DiffractionCorrection\n\n"
                    "\rperturb should be a python list of five floats/ints.\n"
                    "\rYour list contains objects that are not real numbers.\n"
                );
                return -1;
            }
        }

        /*  Store the perturb variable in the DiffractionCorrection class.    */
        tmp = self->perturb;
        Py_INCREF(perturb);
        self->perturb = perturb;
        Py_XDECREF(tmp);
    }

    /*  The input was not a list. Return with error.                          */
    else
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rperturb should be a python list of five floats or ints.\n"
        );
        return -1;
    }

    rssringoccs_Tau_Set_WType(self->wtype, tau);
    rssringoccs_Tau_Set_Psitype(self->psitype, tau);

    /*  Check that all of the following variables are positive using the      *
     *  RSS_RINGOCCSNegValueError preprocessor functions defined above.       */
    if (self->input_res<=0.0)
    {
        RSS_RINGOCCSNegValueError(DiffractionCorrection, res);
        return -1;
    }
    else if (self->sigma<=0.0)
    {
        RSS_RINGOCCSNegValueError(DiffractionCorrection, sigma);
        return -1;
    }
    else if (self->res_factor<=0.0)
    {
        RSS_RINGOCCSNegValueError(DiffractionCorrection, res_factor);
        return -1;
    }

    if (self->ecc<0.0)
    {
        RSS_RINGOCCSNegValueError(DiffractionCorrection, ecc);
        return -1;  
    }
    else
        tau.ecc = self->ecc;

    /*  Check that peri lies in the range -2pi<peri<2pi.                      */
    if (fabs(self->peri)>rssringoccs_Two_Pi)
    {
        RSS_RINGOCCSTwoPiValueError(DiffractionCorrection, peri);
        return -1;
    }
    else
        tau.peri = self->peri;

    /*  Compute the actual resolution to be used for processing.              */
    self->res = self->input_res * self->res_factor;

    /*  If rngreq is still NULL the user didn't supply a range. Set to "all". */
    if (!(rngreq))
        rngreq = PyBytes_FromString("all");

    /*  If the rng variable is a string, make sure it is a legal value and    *
     *  try to extract the corresponding values in kilometers.                */
    if PyBytes_Check(rngreq)
    {
        /*  Convert the Python string to a C string via PyBytes_AsString. */
        char *rng_string = PyBytes_AsString(rngreq);

        /*  Set variable equal to the length of the string.                   */
        unsigned long rng_length = strlen(rng_string);

        /*  Convert the string to lower-case.                                 */
        for (i=0; i<rng_length; ++i)
            rng_string[i] = tolower(rng_string[i]);

        /*  Pass the string into the get range function. This function stores *
         *  the corresponding values in self->range and returns 1 on success, *
         *  and returns 0 on failure. If fails, return with error.            */
        GetRangeFromString(rng_string, self->range);
        if (self->range[0] < 0.0)
        {
            RSS_RINGOCCSRangeError(rng_string)
            return -1;
        }

        /*  Store the rngreq in the diffraction correction class.             */
        tmp = self->rngreq;
        Py_INCREF(rngreq);
        self->rngreq = rngreq;
        Py_XDECREF(tmp);
    }

    /*  If the rng variable is a unicode object (type of string from python)  *
     *  make sure it is a legal value and try to extract the corresponding    *
     *  values in kilometers.                                                 */
    else if PyUnicode_Check(rngreq)
    {
        /*  Convert the Python string to a C string via PyUnicode_AsUTF8. The *
         *  C API recommends not altering the string, so we create a copy of  *
         *  it using strcpy (from string.h).                                  */
        const char *rng_string_cpy = PyUnicode_AsUTF8(rngreq);

        /*  Set variable equal to the length of the string.                   */
        unsigned long rng_length = strlen(rng_string_cpy);

        /*  Allocate enough memory for rng_string so strcpy doesn't fail.     */
        char *rng_string = (char *)malloc(sizeof(char) * rng_length + 1);

        /*  strcpy is a standard library function contained in string.h.      */
        strcpy(rng_string, rng_string_cpy);

        /*  Convert the string to lower-case.                                 */
        for (i=0; i<rng_length; ++i)
            rng_string[i] = tolower(rng_string[i]);

        /*  Pass the string into the get range function. This function stores *
         *  the corresponding values in self->range and returns 1 on success, *
         *  and returns 0 on failure. If fails, return with error.            */
        GetRangeFromString(rng_string, self->range);
        if (self->range[0] < 0.0)
        {
            RSS_RINGOCCSRangeError(rng_string)
            return -1;
        }

        /*  If the above succeeded, we have the range so free rng_string.     */
        else
            free(rng_string);

        /*  Set the rngreq variable in self.                                  */
        tmp = self->rngreq;
        Py_INCREF(rngreq);
        self->rngreq = rngreq;
        Py_XDECREF(tmp);
    }

    /*  If the requested range is a list, try to parse the elements.          */
    else if PyList_Check(rngreq)
    {
        if (!(PyList_Size(rngreq) == RANGE_LEN))
        {
            PyErr_Format(
                PyExc_IndexError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tdiffrec.DiffractionCorrection\n\n"
                "\rrng should be a python list of two floats or ints.\n"
                "\rYour list has %d elements\n", PyList_Size(rngreq)
            );
            return -1;
        }

        iter = PyObject_GetIter(rngreq);
        if (!iter)
        {
            /*  The list is empty, return with an error.                      */
            PyErr_Format(
                PyExc_TypeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tdiffrec.DiffractionCorrection\n\n"
                "\rrng should be a python list of two floats or ints.\n"
            );
            return -1;
        }

        for (i=0; i<RANGE_LEN; ++i)
        {
            next = PyIter_Next(iter);
            if (!next)
            {
                PyErr_Format(
                    PyExc_TypeError,
                    "\n\rError Encountered: rss_ringoccs\n"
                    "\r\tdiffrec.DiffractionCorrection\n\n"
                    "\rrng should be a python list of two floats or ints.\n"
                    "\rSize of your list: %d", i
                );
                return -1;
            }

            /*  Try to parse the elements. Return with error if this fails.   */
            if (PyLong_Check(next))
                self->range[i] = (double)PyLong_AsLong(next);
            else if (PyFloat_Check(next))
                self->range[i] = PyFloat_AsDouble(next);
            else
            {
                PyErr_Format(
                    PyExc_TypeError,
                    "\n\rError Encountered: rss_ringoccs\n"
                    "\r\tdiffrec.DiffractionCorrection\n\n"
                    "\rrng should be a python list of two floats or ints.\n"
                    "\rYour list contains objects that are not real numbers.\n"
                );
                return -1;
            }
        }

        /*  If the user provided values [a, b] with b<a, switch them around.  */
        if (self->range[1] < self->range[0])
        {
            double temp_range_val = self->range[0];
            self->range[0] = self->range[1];
            self->range[1] = temp_range_val;
        }

        /*  Set the rngreq variable.                                          */
        tmp = self->rngreq;
        Py_INCREF(rngreq);
        self->rngreq = rngreq;
        Py_XDECREF(tmp);
    }

    /*  Illegal rng requested. Return with error.                             */
    else
    {
        PyErr_Format(
            PyExc_ValueError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rrng must be a list of two real numbers or a string.\n"
            "\rAllowed strings are:\n"
            "\r\tall               [1.0, 400000.0]\n"
            "\r\tbesselbarnard     [120210.0, 120330.0]\n"
            "\r\tbessel-barnard    [120210.0, 120330.0]\n"
            "\r\tcringripples      [77690.0, 77760.0]\n"
            "\r\tencke             [132900.0, 134200.0]\n"
            "\r\tenckegap          [132900.0, 134200.0]\n"
            "\r\therschel          [118100.0, 118380.0]\n"
            "\r\therschelgap       [118100.0, 118380.0]\n"
            "\r\thuygens           [117650.0, 117950.0]\n"
            "\r\thuygensringlet    [117650.0, 117950.0]\n"
            "\r\tjanusepimetheus   [96200.0, 96800.0]\n"
            "\r\tjeffreys          [118900.0, 119000.0]\n"
            "\r\tjeffreysgap       [118900.0, 119000.0]\n"
            "\r\tkuiper            [119300.0, 119500.0]\n"
            "\r\tkuipergap         [119300.0, 119500.0]\n"
            "\r\tmaxwell           [87410.0, 87610.0]\n"
            "\r\tmaxwellringlet    [87410.0, 87610.0]\n"
            "\r\trussell           [118550.0, 118660.0]\n"
            "\r\trussellgap        [118550.0, 118660.0]\n"
            "\r\ttitan             [77870.0, 77930.0]\n"
            "\r\ttitanringlet      [77870.0, 77930.0\n\n"
        );
        return -1;
    }

    if (self->verbose)
        puts("\tRetrieving variables from DLP instance...");

    /*  Next we're going to run error checks on the input numpy arrays which  *
     *  should be contained inside of the DLPInst object. We'll check that    *
     *  these attributes exist, that they are numpy arrays, are 1 dimensional,*
     *  and have the same number of elements as rho_km_vals. We'll also       *
     *  convert the arrays to double and retrieve a pointer to the data.      *
     *  First, we need to make sure rho_km_vals is a legal numpy array and    *
     *  extract the length of it. Check that rho_km_vals exists in DLPInst.   */
    if (!PyObject_HasAttrString(DLPInst, "rho_km_vals"))
    {
        PyErr_Format(
            PyExc_AttributeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rInput DLP Instance is missing the following attribute:\n"
            "\r\trho_km_vals\n\n"
        );
        return -1;
    }

    /*  If it exists, get a pointer to it.                                    */
    else
        tmp = PyObject_GetAttrString(DLPInst, "rho_km_vals");

    /*  Now make sure rho_km_vals is a numpy array.                           */
    if (!PyArray_Check(tmp))
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rrho_km_vals must be a numpy array.\n"
        );
        return -1;
    }

    /*  If rho_km_vals is a numpy array, try to convert it to double.         */
    else
        arr = PyArray_FromObject(tmp, NPY_DOUBLE, 1, 1);

    /*  If PyArray_FromObject failed arr should be NULL. If so, raise error. */
    if (!arr)
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rCould not convert rho_km_vals to double array. Input is most\n",
            "\rlikely complex numbers or contains a string.\n\n"
        );
        return -1;
    }

    /*  Currently we only allow for one dimensional inputs.                   */
    else if (PyArray_NDIM((PyArrayObject *)arr) != 1)
    {
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rrho_km_vals must be a one-dimensional numpy array.\n"
        );
        return -1;
    }

    /*  If every passed, set tau.rho_km_vals to point to the data inside arr. */
    tau.rho_km_vals = (double *)PyArray_DATA((PyArrayObject *)arr);
    tau.arr_size = PyArray_DIMS((PyArrayObject *)arr)[0];

    /*  Check that there's actually data to process.                          */
    if (tau.arr_size<2)
    {
        PyErr_Format(
            PyExc_ValueError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rrho_km_vals has less than 2 points.\n\n"
        );
        return -1;
    }

    /*  Check all of the variables for errors ensuring they are the correct   *
     *  size and have the right dimension. This is done using the             *
     *  RSS_RINGOCCSCheckArray preprocessor function defined at the start of  *
     *  this file.                                                            */
    RSS_RINGOCCSCheckArray(tau.rho_km_vals, DLPInst, rho_km_vals,
                           tmp, arr, tau.arr_size);
    RSS_RINGOCCSCheckArray(tau.B_rad_vals, DLPInst, B_rad_vals,
                           tmp, arr, tau.arr_size);
    RSS_RINGOCCSCheckArray(tau.phi_rad_vals, DLPInst, phi_rad_vals,
                           tmp, arr, tau.arr_size);
    RSS_RINGOCCSCheckArray(tau.D_km_vals, DLPInst, D_km_vals,
                           tmp, arr, tau.arr_size);
    RSS_RINGOCCSCheckArray(tau.f_sky_hz_vals, DLPInst, f_sky_hz_vals,
                           tmp, arr, tau.arr_size);
    RSS_RINGOCCSCheckArray(tau.rho_dot_kms_vals, DLPInst,
                           rho_dot_kms_vals, tmp, arr, tau.arr_size);
    RSS_RINGOCCSCheckArray(tau.t_oet_spm_vals, DLPInst,
                           t_oet_spm_vals, tmp, arr, tau.arr_size);
    RSS_RINGOCCSCheckArray(tau.t_ret_spm_vals, DLPInst,
                           t_ret_spm_vals, tmp, arr, tau.arr_size);
    RSS_RINGOCCSCheckArray(tau.t_set_spm_vals, DLPInst,
                           t_set_spm_vals, tmp, arr, tau.arr_size);
    RSS_RINGOCCSCheckArray(tau.rho_corr_pole_km_vals, DLPInst,
                           rho_corr_pole_km_vals,tmp, arr, tau.arr_size);
    RSS_RINGOCCSCheckArray(tau.rho_corr_timing_km_vals, DLPInst,
                           rho_corr_timing_km_vals, tmp, arr, tau.arr_size);
    RSS_RINGOCCSCheckArray(tau.raw_tau_threshold_vals, DLPInst,
                           raw_tau_threshold_vals, tmp, arr, tau.arr_size);
    RSS_RINGOCCSCheckArray(tau.phi_rl_rad_vals, DLPInst,
                           phi_rl_rad_vals, tmp, arr, tau.arr_size);
    RSS_RINGOCCSCheckArray(tau.p_norm_vals, DLPInst, p_norm_vals,
                           tmp, arr, tau.arr_size);
    RSS_RINGOCCSCheckArray(tau.phase_rad_vals, DLPInst,
                           phase_rad_vals, tmp, arr, tau.arr_size);

    if (!PyObject_HasAttrString(DLPInst, "history"))
    {
        PyErr_Format(
            PyExc_AttributeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rInput DLP Instance is missing the following attribute:\n"
            "\r\thistory\n\n"
        );
        return -1;   
    }
    else
        history = PyObject_GetAttrString(DLPInst, "history");

    /*  Store the dlp history inside of the DiffractionCorrection class. This *
     *  tmp, Py_INCREF, Py_XDECREF method is recommended in the Python C-API  *
     *  documentation as a means of safely storing the variable.              */
    tmp = self->dathist;
    Py_INCREF(history);
    self->dathist = history;
    Py_XDECREF(tmp);

    RSS_RINGOCCSCheckArrayPositive(tau.D_km_vals, tau.arr_size);
    RSS_RINGOCCSCheckArrayPositive(tau.rho_km_vals, tau.arr_size);
    RSS_RINGOCCSCheckArrayPositive(tau.f_sky_hz_vals, tau.arr_size);
    RSS_RINGOCCSCheckArrayTwoPi(tau.B_rad_vals, tau.arr_size);
    RSS_RINGOCCSCheckArrayTwoPi(tau.phi_rad_vals, tau.arr_size);

    /*  Make sure the power isn't negative.                                   */
    if (rssringoccs_Min_Double(tau.p_norm_vals, tau.arr_size) < 0.0)
    {
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\ttau.p_norm_vals contains non-positive values.\n"
        );
        return -1;
    }

    /*  Run an error check as to whether or not this is an egress, ingress,   *
     *  or chord occultation, and if the data can be processed accordingly.   */
    double min_val = rssringoccs_Min_Double(tau.rho_dot_kms_vals, tau.arr_size);
    double max_val = rssringoccs_Max_Double(tau.rho_dot_kms_vals, tau.arr_size);
    double dx_km   = tau.rho_km_vals[1] - tau.rho_km_vals[0];

    if (dx_km == 0.0)
    {
        PyErr_Format(
            PyExc_ValueError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rdx_km is zero: rho_km_vals[1]-rho_km_vals[0] = 0.\n\n"
        );
        return -1;
    }
    else if (self->res < 1.99*dx_km)
    {
        PyErr_Format(
            PyExc_ValueError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rResolution is less than twice the sample space.\n"
            "\r\tRequested resolution (km):     %f\n"
            "\r\tSample spacing (km):           %f\n\n"
            "\r\tChoose a resolution GREATER than %f (km)\n\n"
            "\rPLEASE NOTE:\n"
            "\r\tTo be consistent with PDS results, a factor of 0.75 is\n"
            "\r\tapplied to the requested resolution. To ignore set, set\n"
            "\r\tthe 'res_factor=0' when calling DiffractionCorrection.\n"
            "\r\t\tres_factor currently set to:    %f",
            self->res, dx_km, 2.0*dx_km/self->res_factor, self->res_factor
        );
        return -1;
    }

    if ((min_val < 0.0) && (max_val > 0.0))
    {
        PyErr_Format(
            PyExc_ValueError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r\tdrho/dt has positive and negative values.\n"
            "\r\tYour input file is probably a chord occultation.\n"
            "\r\tDiffraction Correction can only be performed for\n"
            "\r\tone event at a time. That is, ingress or egress.\n\n"
            "\r\tTO CORRECT THIS:\n"
            "\r\t\tSplit the input into two parts: Ingress and Engress\n"
            "\r\t\tand perform diffraction correction twice.\n\n"
        );
        return -1;
    }
    else if ((min_val == 0.0) || (max_val == 0.0))
    {
        PyErr_Format(
            PyExc_ValueError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r\tdrho/dt has zero valued elements.\n"
            "\r\tYour input file is probably a chord occultation.\n"
            "\r\tDiffraction Correction can only be performed for\n"
            "\r\tone event at a time. That is, ingress or egress.\n\n"
            "\r\tTO CORRECT THIS:\n"
            "\r\t\tSplit the input into two parts: Ingress and Engress\n"
            "\r\t\tand perform diffraction correction twice.\n\n"
        );
        return -1;
    }
    else if ((dx_km>0.0) && (max_val<0.0))
    {
        for(i=0; i<tau.arr_size; ++i)
            tau.rho_dot_kms_vals[i] = fabs(tau.rho_dot_kms_vals[i]);
    }
    else if ((dx_km < 0.0) && (min_val > 0.0))
    {
        PyErr_Format(
            PyExc_ValueError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r\trho_km_vals is decreasing yet rho_dot_kms_vals\n"
            "\r\tis positiive. Check DLP class for errors.\n\n"
        );
        return -1;
    }
    else if (dx_km<0.0)
    {
        rssringoccs_Reverse_Double_Array(tau.rho_km_vals,      tau.arr_size);
        rssringoccs_Reverse_Double_Array(tau.phi_rad_vals,     tau.arr_size);
        rssringoccs_Reverse_Double_Array(tau.B_rad_vals,       tau.arr_size);
        rssringoccs_Reverse_Double_Array(tau.D_km_vals,        tau.arr_size);
        rssringoccs_Reverse_Double_Array(tau.phase_rad_vals,   tau.arr_size);
        rssringoccs_Reverse_Double_Array(tau.p_norm_vals,      tau.arr_size);
        rssringoccs_Reverse_Double_Array(tau.f_sky_hz_vals,    tau.arr_size);
        rssringoccs_Reverse_Double_Array(tau.rho_dot_kms_vals, tau.arr_size);
        for(i=0; i<tau.arr_size; ++i)
            tau.rho_dot_kms_vals[i] = fabs(tau.rho_dot_kms_vals[i]);
        dx_km *= -1.0;
    }

    if (self->verbose)
        puts("\tComputing Necessary Variables...");

    double *lambda_sky = (double *)malloc(sizeof(double) * tau.arr_size);
    double *mu_vals    = (double *)malloc(sizeof(double) * tau.arr_size);

    /*  Allocate memory for the diffracted data and the fresnel scale.        */
    tau.T_in      = malloc(sizeof(*tau.T_in)      * tau.arr_size);
    tau.F_km_vals = malloc(sizeof(*tau.F_km_vals) * tau.arr_size);

    for (i=0; i<tau.arr_size; ++i)
    {
        /*  Compute the complex amplitude, T_hat_vals.                        */
        tau.T_in[i] = rssringoccs_CDouble_Polar(sqrt(tau.p_norm_vals[i]),
                                                tau.phase_rad_vals[i]);

        /*  Compute the wavelength lambda and scale-factor mu.                */
        lambda_sky[i] = rssringoccs_Double_Frequency_To_Wavelength(tau.f_sky_hz_vals[i]);

        mu_vals[i]    = sin(fabs(tau.B_rad_vals[i]));

        /*  And finally, compute the Fresnel scale.                           */
        tau.F_km_vals[i]  = rssringoccs_Double_Fresnel_Scale(lambda_sky[i],
                                                             tau.D_km_vals[i],
                                                             tau.phi_rad_vals[i],
                                                             tau.B_rad_vals[i]);
    }

    /*  Get the normalized equivalent width of the window type.               */
    GetNormeqFromString(tau.wtype, &self->norm_eq);

    /*  Use calloc to both allocate memory for tau.w_km_vals (like malloc)    *
     *  and initialize the data to zero (unlike malloc). This is similar to   *
     *  numpy.zeros(tau.arr_size) in Python.                                  */
    tau.w_km_vals = calloc(tau.arr_size, sizeof(double));

    /*  Declare long pointer-to-pointer which stores the indices where        *
     *  F_km_vals is non-zero in the first slot (Prange[0]), and the size of  *
     *  this array in the second (*Prange[1]).                                */
    unsigned long **Prange;
    unsigned long *Prange_Index;
    unsigned long Prange_Size;
    double w_fac;

    if (self->bfac)
    {
        w_fac = self->norm_eq;
        double omega;
        double *alpha  = malloc(sizeof(*alpha) * tau.arr_size);
        double *P_vals = malloc(sizeof(*P_vals) * tau.arr_size);

        for(i=0; i<tau.arr_size; ++i)
        {
            omega      = rssringoccs_Two_Pi * tau.f_sky_hz_vals[i];
            alpha[i]   = omega * self->sigma;
            alpha[i]  *= alpha[i] * 0.5 / tau.rho_dot_kms_vals[i];
            P_vals[i]  = self->res/(alpha[i]*tau.F_km_vals[i]*tau.F_km_vals[i]);
        }

        Prange = rssringoccs_Where_Greater_Double(P_vals, tau.arr_size, 1.0);
        Prange_Index = Prange[0];
        Prange_Size  = *Prange[1];

        for (i=0; i<Prange_Size; ++i)
            tau.w_km_vals[Prange_Index[i]] = w_fac *
                rssringoccs_Double_Resolution_Inverse(P_vals[Prange_Index[i]]) / alpha[i];

        free(P_vals);
        free(alpha);
    }
    else
    {
        w_fac = self->norm_eq/self->res;

        for (i=0; i<tau.arr_size; ++i)
            tau.w_km_vals[i] = 2.0*tau.F_km_vals[i]*tau.F_km_vals[i]*w_fac;

        Prange = rssringoccs_Where_Greater_Double(tau.F_km_vals, tau.arr_size, 0.0);
        Prange_Index = Prange[0];
        Prange_Size  = *Prange[1];
    }

    double *rho_legal = malloc(sizeof(*rho_legal) * Prange_Size);

    for(i=0; i<Prange_Size; ++i)
        rho_legal[i] = tau.rho_km_vals[Prange_Index[i]] -
                       0.5*tau.w_km_vals[Prange_Index[i]];

    max_val = rssringoccs_Max_Double(rho_legal, Prange_Size);

    for(i=0; i<Prange_Size; ++i)
        rho_legal[i] = tau.rho_km_vals[Prange_Index[i]] +
                       0.5*tau.w_km_vals[Prange_Index[i]];
    min_val = rssringoccs_Min_Double(rho_legal, Prange_Size);

    unsigned long **wrange = rssringoccs_Where_LesserGreater_Double(tau.rho_km_vals, tau.arr_size,
                                                                    min_val, max_val);
    unsigned long *wrange_Index = wrange[0];
    unsigned long wrange_Size = *wrange[1];


    if (wrange_Size == 0)
    {
        PyErr_Format(
            PyExc_ValueError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r\tThe window width is too large to reconstruct anything.\n"
        );
        return -1;
    }
    else if (max_val < self->range[0])
    {
        PyErr_Format(
            PyExc_ValueError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r\tMinimum requested range is greater than available data.\n"
            "\r\t\tYour Requested Minimum (km):    %lld\n"
            "\r\t\tYour Requested Maximum (km):    %lld\n"
            "\r\t\tMaximum Available Data (km):    %lld\n",
            (long long)self->range[0], (long long)self->range[1],
            (long long)max_val
        );
        return -1;
    }
    else if (min_val > self->range[1])
    {
        PyErr_Format(
            PyExc_ValueError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r\tMinimum requested range is greater than available data.\n"
            "\r\t\tYour Requested Minimum:  %lld km\n"
            "\r\t\tYour Requested Maximum:  %lld km\n"
            "\r\t\tMinimum Available Data:  %lld km\n",
            (long long)self->range[0], (long long)self->range[1],
            (long long)max_val
        );
        return -1;
    }
    else
    {
        free(wrange[0]);
        free(wrange[1]);
        free(wrange);

        wrange = rssringoccs_Where_LesserGreater_Double(tau.rho_km_vals, tau.arr_size,
                                                        min_val, max_val);
        wrange_Index = wrange[0];
        wrange_Size = *wrange[1];

        self->start = wrange_Index[0];
        self->finish = wrange_Index[wrange_Size-1];
        self->n_used = self->finish - self->start;
    }

    tau.start  = self->start;
    tau.n_used = self->n_used;

    /*  For the first inverse calculatiion, tau.use_fwd must be set to false. */
    tau.use_fwd = self->use_fft;
    tau.use_norm = self->use_norm;

    tau.kd_vals = (double *)malloc(sizeof(double) * tau.arr_size);
    for(i=0; i<tau.arr_size; ++i)
        tau.kd_vals[i] = rssringoccs_Double_Wavelength_To_Wavenumber(lambda_sky[i]) *
                     tau.D_km_vals[i];

    /*  Free the variables you don't need anymore.                            */
    free(wrange_Index);
    free(wrange);
    free(Prange_Index);
    free(Prange);
    free(rho_legal);
    free(lambda_sky);

    if (tau.start > tau.arr_size){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rStarting index (start) is greater than the size of the array.\n"
        );
        return -1;
    }
    else if (tau.start+tau.n_used > tau.arr_size){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rFinal index (start+n_used) is greater than size of array.\n"
        );
        return -1;
    }

    tau.T_out = malloc((tau.n_used+1)*sizeof(*tau.T_out));
    /*
    if (self->use_fft)
        DiffractionCorrectionSimpleFFT(&tau);
    else
    {
        if (tau.order == 0)
        {
            if ((tau.ecc == 0.0) && (tau.peri == 0.0))
                if ((tau.perturb[0] == 0) && (tau.perturb[1] == 0) &&
                    (tau.perturb[2] == 0) && (tau.perturb[3] == 0) &&
                    (tau.perturb[4] == 0))
                    DiffractionCorrectionNewton(&tau);
                else
                    DiffractionCorrectionPerturbedNewton(&tau);
            else
                DiffractionCorrectionEllipse(&tau);
        }
        else if (tau.order == 1)
            DiffractionCorrectionFresnel(&tau);
        else
            DiffractionCorrectionLegendre(&tau);
    }
    */
    time_t t1, t2;

    t1 = clock();
    if (strcmp(tau.psitype, "fresnel") == 0)
        DiffractionCorrectionFresnel(&tau);
    else if (strcmp(tau.psitype, "newton") == 0)
    {
        if ((tau.perturb[0] == 0) && (tau.perturb[1] == 0) &&
            (tau.perturb[2] == 0) && (tau.perturb[3] == 0) &&
            (tau.perturb[4] == 0))
            DiffractionCorrectionNewton(&tau);
        else
            DiffractionCorrectionPerturbedNewton(&tau);

    }
    else
        DiffractionCorrectionNewton(&tau);
    t2 = clock();

    printf("%f\n", (double)(t2-t1) / CLOCKS_PER_SEC);

    if (tau.error_occurred)
    {
        PyErr_Format(PyExc_RuntimeError, tau.error_message);
        return -1;
    }
    else
        tau.n_used += 1;

    RSS_RINGOCCSSetArray(self, tau, rho_km_vals, tmp, NPY_DOUBLE,
                         tau.start, tau.n_used);

    long pynused = (long)tau.n_used;
    PyObject *T_out = PyArray_SimpleNewFromData(1, &pynused,
                                                NPY_CDOUBLE, (void *)tau.T_out);

    tmp = self->T_vals;
    Py_INCREF(T_out);
    self->T_vals = T_out;
    Py_XDECREF(tmp);

    PyObject *T_in = PyArray_SimpleNewFromData(1, &pynused,
                                               NPY_CDOUBLE, (void *)tau.T_in);

    tmp = self->T_hat_vals;
    Py_INCREF(T_in);
    self->T_hat_vals = T_in;
    Py_XDECREF(tmp);

    /*  Reset the phase to normal. This is very important! We modified the    *
     *  data in the dlp_inst directly, rather than steal a copy, so we need   *
     *  to put it back. Luckily, we only negated it. So just re-negate it.    */
    for (i=0; i<tau.arr_size; ++i)
        tau.phase_rad_vals[i] = -tau.phase_rad_vals[i];

    /*  New free everything we've malloc'd, or we'll have memory leaks. Do    *
     *  NOT free the stuff malloc'd that we are storing in self. This will be *
     *  freed once the instance of the DiffractionCorrection class is         *
     *  destroyed. Only free the C-level stuff.                               */
    free(tau.wtype);
    free(tau.psitype);
    free(tau.kd_vals);

    return 1;
}

static PyObject *
DiffractionCorrection(PyDiffrecObj *self, PyObject *Py_UNUSED(ignored))
{
    return PyUnicode_FromFormat("DiffractionCorrection");
}

static PyMethodDef DiffractionCorrection_methods[] =
{
    {
        "DiffractionCorrection",
        (PyCFunction) DiffractionCorrection,
        METH_NOARGS,
        "Diffraction correction class."
    },
    {NULL}
};

static PyMemberDef Custom_members[] = {
    {"rho_km_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, rho_km_vals), 0,
     "Ring radius"},
    {"phase_rad_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, phase_rad_vals), 0,
     "Phase"},
    {"B_rad_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, B_rad_vals), 0,
     "Ring inclination angle"},
    {"D_km_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, D_km_vals), 0,
     "Spacecraft to ring-intercept point distance"},
    {"F_km_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, F_km_vals), 0,
     "Fresnel scale"},
    {"T_hat_fwd_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, T_hat_fwd_vals), 0,
     "Forward modeling of data"},
    {"T_hat_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, T_hat_vals), 0,
     "Raw diffraction data"},
    {"T_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, T_vals), 0,
     "Diffraction corrected data"},
    {"f_sky_hz_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, f_sky_hz_vals), 0,
     "Frequency of the input signal"},
    {"outfiles", T_OBJECT_EX, offsetof(PyDiffrecObj, outfiles), 0,
     "CSV of corrected data"},
    {"p_norm_fwd_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, p_norm_fwd_vals), 0,
     "Forward modeling of power"},
    {"p_norm_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, p_norm_vals), 0,
     "Raw power data"},
    {"phase_fwd_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, phase_fwd_vals), 0,
     "Forward modeling of phase"},
    {"phase_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, phase_vals), 0,
     "Diffraction corrected phase"},
    {"phi_rad_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, phi_rad_vals), 0,
     "Ring azimuth angle"},
    {"phi_rl_rad_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, phi_rl_rad_vals), 0,
     "Ring longitude angle"},
    {"power_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, power_vals), 0,
     "Diffraction corrected power"},
    {"raw_tau_threshold_vals", T_OBJECT_EX,
     offsetof(PyDiffrecObj, raw_tau_threshold_vals), 0,
     "Raw threshold optical depth"},
    {"rev_info", T_OBJECT_EX, offsetof(PyDiffrecObj, rev_info), 0,
     "Information about the occultation"},
    {"rho_corr_pole_km_vals", T_OBJECT_EX,
     offsetof(PyDiffrecObj, rho_corr_pole_km_vals), 0,
     "Ring radius with pole correction."},
    {"rho_corr_timing_km_vals", T_OBJECT_EX,
     offsetof(PyDiffrecObj, rho_corr_timing_km_vals), 0,
     "Ring radius with timing correction."},
    {"rho_dot_kms_vals", T_OBJECT_EX,
     offsetof(PyDiffrecObj, rho_dot_kms_vals), 0,
     "Time derivative of the ring radius."},
    {"t_oet_spm_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, t_oet_spm_vals), 0,
     "Observed event time in seconds past midnight"},
    {"t_ret_spm_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, t_ret_spm_vals), 0,
     "Ring event time in seconds past midnight"},
    {"t_set_spm_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, t_set_spm_vals), 0,
     "Spacecraft event time in seconds past midnight"},
    {"tau_threshold_vals", T_OBJECT_EX,
     offsetof(PyDiffrecObj, tau_threshold_vals), 0,
     "Diffraction corrected threshold optical depth"},
    {"tau_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, tau_vals), 0,
     "Optical depth"},
    {"w_km_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, w_km_vals), 0,
     "window width as a function of ring radius"},
    {"dathist", T_OBJECT_EX, offsetof(PyDiffrecObj, dathist), 0,
     "History of input dlp instance"},
    {"history", T_OBJECT_EX, offsetof(PyDiffrecObj, history), 0,
     "History of this tau instance"},
    {"perturb", T_OBJECT_EX, offsetof(PyDiffrecObj, perturb), 0,
     "Coefficients for perturbation polynomial"},
    {"rngreq", T_OBJECT_EX, offsetof(PyDiffrecObj, rngreq), 0,
     "Range requested by user"},
    {"use_fft", T_BOOL, offsetof(PyDiffrecObj, use_fft), 0,
     "Use of FFTs for data processing"},
    {"bfac", T_BOOL, offsetof(PyDiffrecObj, bfac), 0,
     "Use of b-factor in window width"},
    {"verbose", T_BOOL, offsetof(PyDiffrecObj, verbose), 0,
     "Print status updates"},
    {"use_norm", T_BOOL, offsetof(PyDiffrecObj, use_norm), 0,
     "Use of window normalization"},
    {"fwd", T_BOOL, offsetof(PyDiffrecObj, fwd), 0,
     "Forward modeling Boolean"},

    {NULL}  /* Sentinel */
};

static PyTypeObject DiffrecType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "diffrec.DiffractionCorrection",
    .tp_doc = "Diffraction Correction class.",
    .tp_basicsize = sizeof(PyDiffrecObj),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = PyType_GenericNew,
    .tp_init = (initproc) Diffrec_init,
    .tp_dealloc = (destructor) Diffrec_dealloc,
    .tp_members = Custom_members,
    .tp_methods = DiffractionCorrection_methods,
};

static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    .m_name = "custom",
    .m_doc = "Example module that creates an extension type.",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit_diffrec(void)
{
    PyObject *m;
    if (PyType_Ready(&DiffrecType) < 0)
        return NULL;

    m = PyModule_Create(&moduledef);

    if (m == NULL)
        return NULL;

    Py_INCREF(&DiffrecType);
    if (PyModule_AddObject(m, "DiffractionCorrection",
                           (PyObject *) &DiffrecType) < 0)
    {
        Py_DECREF(&DiffrecType);
        Py_DECREF(m);
        return NULL;
    }

    import_array();
    return m;
}
