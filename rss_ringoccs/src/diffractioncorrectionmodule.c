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
 *  ReverseDoubleArray:                                                       *
 *      The standard O(n) algorithm for taking an array and reversing its     *
 *      order. This is equivalent to arr = arr[::-1] for a numpy array arr    *
 *      in Python 3.                                                          *
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

/*  Support for booleans was added in C99, and various definitions are found  *
 *  in stdbool.h.                                                             */
#include <stdbool.h>

/*  Support for complex variables was added in C99. This file contains the    *
 *  definition of complex data types and convenient functions for             *
 *  manipulating the data.                                                    */
#include <complex.h>

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
#include <rss_ringoccs.h>
#include <rss_ringoccs_math_constants.h>
#include <rss_ringoccs_special_functions.h>

/*  This header file contains pure C99 code, but relies on the NON-STANDARD   *
 *  FFTW3 library for performing FFTs. As part of the installation, you will  *
 *  need to install FFTW. The setup scripts will link this library            *
 *  automatically without user interaction. This contains all of the          *
 *  diffraction correction routines which return the corrected data.          */
#include "_diffraction_correction.h"
#include "_diffraction_modeling.h"

/*  All of the "psi" functions, included Newton-Raphson routines, are found   *
 *  here. It uses C99 only, and no non-standard libraries are needed.         */
#include "_fresnel_kernel.h"


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
    long i;                                                                    \
    outtype *out = (outtype *)malloc(sizeof(outtype)*dim);                     \
    for (i=0; i<dim; ++i)                                                      \
        out[i] = (*f)(((intype *)x1)[i], x2, x3);                              \
    y = out;                                                                   \
})

#define _array_from_four(x1, x2, x3, x4, y, dim, intype, outtype, f) ({        \
    long i;                                                                    \
    outtype *out = (outtype *)malloc(sizeof(outtype)*dim);                     \
    for (i=0; i<dim; ++i)                                                      \
        out[i] = (*f)(((intype *)x1)[i], x2, x3, x4);                          \
    y = out;                                                                   \
})

/*  The Python wrappers for several of the functions is completely identical  *
 *  with the exception of which C function needs to be called, and what the   *
 *  function will be called at the Python level. We create a macro with two   *
 *  inputs:                                                                   *
 *      FuncName:                                                             *
 *          The name of the function at the python level.                     *
 *      CName:                                                                *
 *          The name of the function at the C level with the data type        *
 *          excluded. For example, if you wanted to wrap the Sinc_Double      *
 *          function defined in special_functions/, CName whould be Sinc.     *
 *  NOTE:                                                                     *
 *      Functions built using this macro MUST require four variables of the   *
 *      following types:                                                      *
 *          rho:                                                              *
 *              A numpy array or list of real numbers, a float, or an int.    *
 *          a:                                                                *
 *              A positive real number.                                       *
 *          b:                                                                *
 *              A positive real number GREATER than a.                        *
 *          F:                                                                *
 *              A positive real number.                                       *
 *      Ideally, only use this for the diffraction modeling functions which   *
 *      require four inputs. rho should be ring radius, a and b should be     *
 *      specific radii, and F is the fresnel scale. rho, a, b, and F should   *
 *      all have the same units, preferably kilometers (km).                  */
#define RSS_RINGOCCSDiffModelingFourVars(FuncName, CName)                      \
static PyObject *FuncName(PyObject *self, PyObject *args)                      \
{                                                                              \
    /*  We'll need output and capsule for safely creating the output array   */\
    /*  and ensuring we don't have a memory leak. rho is the input variable. */\
    PyObject *output, *capsule,  *rho;                                         \
                                                                               \
    /*  a, b, and F are the input real numbers. Currently the code only      */\
    /*  accepts floats and ints, and does not deal with fancy Python objects.*/\
    double a, b, F;                                                            \
                                                                               \
    /*  Variable for the size of the input array which we'll retrieve later. */\
    long dim;                                                                  \
                                                                               \
    /*  Try to extract the values passed from the user, return error if fail.*/\
    if (!PyArg_ParseTuple(args, "Oddd", &rho, &a, &b, &F))                     \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rCould not parse inputs. Legal inputs are:\n"                    \
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"       \
            "\r\ta:     Positive constant (Float)\n"                           \
            "\r\tb:     Positive constant (Float) greater than a\n"            \
            "\r\tF      Positive constant (Float)\n\n"                         \
            "\rNotes:\n"                                                       \
            "\r\trho must be a non-empty one dimensional numpy array,\n"       \
            "\r\ta float, or an int.\n", VarToString(FuncName)                 \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  For all functions using this macro, a <= b. Check this.              */\
    if (a >= b)                                                                \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_ValueError,                                                  \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInner radius is not less than outer radius (i.e. a >= b).\n",   \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  We're modeling in polar coordinates, so radius values should not be  */\
    /*  negative. Check this and raise an error if a < 0.                    */\
    else if (a < 0.0)                                                          \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInner radius is not positive. (i.e. a<0)\n",                    \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  The Fresnel scale should never be zero or negative. Check.           */\
    else if (F <= 0.0)                                                         \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rFresnel scale is not positive (i.e. F<=0).\n",                  \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  If the user supplied a float or int, simply pass the inputs to the   */\
    /*  C function and return the value back to Python.                      */\
    if (PyLong_Check(rho) || PyFloat_Check(rho))                               \
    {                                                                          \
        /*  Convert a python float or int into a C double.                   */\
        double val = PyFloat_AsDouble(rho);                                    \
                                                                               \
        /*  The function returns a complex variable, so declare one.         */\
        complex double out;                                                    \
                                                                               \
        /*  Now get the computed value and convert the C complex double to a */\
        /*  Python complex value and return. The ## symbol tells the         */\
        /*  preprocessor to concatenate the words for the compiler. So if we */\
        /*  pass Gap_Diffraction as CName, the compiler will see             */\
        /*      out = Gap_Diffraction_Double(val, a, b, F);                  */\
        /*  which is exactly what we want.                                   */\
        out = CName##_Double(val, a, b, F);                                    \
                                                                               \
        /*  To build a Python complex value we need to pass the real and     */\
        /*  imaginary parts of out into PyComplex_FromDoubles. We can easily */\
        /*  access these values using creal and cimag, which are functions   */\
        /*  found in complex.h.                                              */\
        return PyComplex_FromDoubles(creal(out), cimag(out));                  \
    }                                                                          \
                                                                               \
    /*  As of v1.3 we allow lists to be passed to these routines.            */\
    else if (PyList_Check(rho))                                                \
    {                                                                          \
        /*  If the user passed a list, we'll need to return one. Create a    */\
        /*  variable for indexing over the input list.                       */\
        long i;                                                                \
                                                                               \
        /*  Create another variable for indexing over the list. This will be */\
        /*  the object corresponding to the value of the index i.            */\
        PyObject *ith_item;                                                    \
                                                                               \
        /*  And lastly, a double for storing the input value and a complex   */\
        /*  double for storing the output value.                             */\
        double val;                                                            \
        double complex out;                                                    \
                                                                               \
        /*  Get the number of elements in the list.                          */\
        dim = PyList_Size(rho);                                                \
                                                                               \
        /*  If the input list is empty, return with an error.                */\
        if (dim == 0)                                                          \
        {                                                                      \
            PyErr_Format(                                                      \
                PyExc_TypeError,                                               \
                "\n\rError Encountered: rss_ringoccs\n"                        \
                "\r\tdiffrec.%s\n\n"                                           \
                "\rInput list is empty.\n", VarToString(FuncName)              \
            );                                                                 \
            return NULL;                                                       \
        }                                                                      \
                                                                               \
        /*  Create a new list to be returned to the user.                    */\
        output = PyList_New(dim);                                              \
                                                                               \
        /*  Loop over the entries of the list, examine values for error, and */\
        /*  try to perform computations if legal inputs are given.           */\
        for (i=0; i<dim; ++i)                                                  \
        {                                                                      \
            /*  Extract the python object in the ith slot of the list.       */\
            ith_item = PyList_GET_ITEM(rho, i);                                \
                                                                               \
            /*  The list should be homogeneous with real numbers only. Check */\
            /*  each item and raise an error otherwise.                      */\
            if (!PyFloat_Check(ith_item) && !PyLong_Check(ith_item))           \
            {                                                                  \
                PyErr_Format(PyExc_TypeError,                                  \
                             "\n\rError Encountered: rss_ringoccs\n"           \
                             "\r\tdiffrec.%s\n\n"                              \
                             "\rInput list must contain real numbers only.\n", \
                             VarToString(FuncName));                           \
                return NULL;                                                   \
            }                                                                  \
                                                                               \
            /*  Convert the python float to a C double.                      */\
            val = PyFloat_AsDouble(ith_item);                                  \
                                                                               \
            /*  Compute the value, set the item in the list, and move on to  */\
            /*  to the next one. As explained before, ## concatenates words  */\
            /*  so CName##_Double compiles as CName_Double.                  */\
            out = CName##_Double(val, a, b, F);                                \
            PyList_SET_ITEM(output, i,                                         \
                            PyComplex_FromDoubles(creal(out), cimag(out)));    \
        }                                                                      \
                                                                               \
        /*  The list created without errors, so just return it.              */\
        return output;                                                         \
    }                                                                          \
                                                                               \
    /*  If the input was not a list or a number, it must be a numpy array.   */\
    /*  Check this and return error otherwise.                               */\
    else if (!PyArray_Check(rho))                                              \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rCould not parse inputs. Legal inputs are:\n"                    \
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"       \
            "\r\ta:     Positive constant (Float)\n"                           \
            "\r\tb:     Positive constant (Float) greater than a\n"            \
            "\r\tF      Positive constant (Float)\n\n"                         \
            "\rNotes:\n"                                                       \
            "\r\trho must be a non-empty one dimensional numpy array,\n"       \
            "\r\ta float, or an int.\n", VarToString(FuncName)                 \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  Check the input array to make sure it is valid. It is important to   */\
    /*  check the number of dimensions BEFORE setting the dim variable. If   */\
    /*  the input array is empty, PyArray_DIMS(rho) will be a NULL pointer.  */\
    /*  Trying to access it with dim = PyArray_DIMS(rho)[0] will create a    */\
    /*  segmentation fault, crashing the Python interpreter.                 */\
    if (PyArray_NDIM((PyArrayObject *)rho) != 1)                               \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInput numpy array is not one-dimensional.\n",                   \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  Useful information about the data. The numpy API functions want a    */\
    /*  PyArrayObject pointer, so we need to cast rho as this type.          */\
    int typenum = PyArray_TYPE((PyArrayObject *)rho);                          \
    void *data  = PyArray_DATA((PyArrayObject *)rho);                          \
    dim         = PyArray_DIMS((PyArrayObject *)rho)[0];                       \
                                                                               \
    /*  This variable is for the data type of the output variable. It should */\
    /*  the complex valued equivalent of the input. So NPY_FLOAT should      */\
    /*  return NPY_CLOAT. All integer types will return NPY_CDOUBLE.         */\
    /*  If the input array was an int-like array (long, short, char, etc)    */\
    /*  then the output will be complex double by default/                   */\
    int outtype;                                                               \
                                                                               \
    /*  This void pointer will point to the output data we'll create later.  */\
    void *y;                                                                   \
                                                                               \
    /*  If the input array is empty, return with an error.                   */\
    if (dim == 0)                                                              \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInput numpy array is empty.\n", VarToString(FuncName)           \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  For float, double, and long double precision (which corresponds to   */\
    /*  numpy.float32, numpy.float64, and numpy.float128, respectively) we   */\
    /*  have functions which compute at these respective precisions. Pass    */\
    /*  the values into the _array_from_four macro which will loop over the  */\
    /*  data and create a new output pointer, which y will then point to.    */\
    if (typenum == NPY_FLOAT)                                                  \
    {                                                                          \
        _array_from_four(data, a, b, F, y, dim, float,                         \
                         complex float, CName##_Float);                        \
        outtype = NPY_CFLOAT;                                                  \
    }                                                                          \
    else if (typenum == NPY_DOUBLE)                                            \
    {                                                                          \
        _array_from_four(data, a, b, F, y, dim, double,                        \
                         complex double, CName##_Double);                      \
        outtype = NPY_CDOUBLE;                                                 \
    }                                                                          \
    else if (typenum == NPY_LONGDOUBLE)                                        \
    {                                                                          \
        _array_from_four(data, a, b, F, y, dim, long double,                   \
                         complex long double, CName##_Long_Double);            \
        outtype = NPY_CLONGDOUBLE;                                             \
    }                                                                          \
                                                                               \
    /*  For all other types we try to convert to double and compute.         */\
    else                                                                       \
    {                                                                          \
                                                                               \
        /*  Try to convert the input numpy array to double and compute.      */\
        PyObject *newrho = PyArray_FromObject(rho, NPY_DOUBLE, 1, 1);          \
                                                                               \
        /*  If PyArray_FromObject failed, newrho should be NULL. Check this. */\
        if (!(newrho))                                                         \
        {                                                                      \
            PyErr_Format(PyExc_TypeError,                                      \
                        "\n\rError Encountered: rss_ringoccs\n"                \
                        "\r\tdiffrec.%s\n\n"                                   \
                        "\rInvalid data type for input array. Input should"    \
                        "\n\rbe a 1-dimensional array of real numbers.\n",     \
                        VarToString(FuncName));                                \
            return NULL;                                                       \
        }                                                                      \
        else                                                                   \
            _array_from_four(data, a, b, F, y, dim, double,                    \
                             complex double, CName##_Double);                  \
                                                                               \
        /*  Retrieve a pointer to the data inside of newrho and set outtype. */\
        data = PyArray_DATA((PyArrayObject *)newrho);                          \
        outtype = NPY_CDOUBLE;                                                 \
    }                                                                          \
                                                                               \
    /*  Create a numpy array from the newly computed data and set a capsule  */\
    /*  for it. This ensures when the corresponding Python variable is       */\
    /*  deleted or removed the memory allocated to the respective C pointer  */\
    /*  is freed. Skipping this will result in memory leaks!                 */\
    output  = PyArray_SimpleNewFromData(1, &dim, outtype, y);                  \
    capsule = PyCapsule_New(y, NULL, capsule_cleanup);                         \
                                                                               \
    /*  This frees the variable at the Python level once it's destroyed.     */\
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);                   \
                                                                               \
    /*  Return the results to Python.                                     */   \
    return Py_BuildValue("N", output);                                         \
}

/*  Same as RSS_RINGOCCSDiffModelingFourVars but for when three variables are *
 *  needed instead of four.                                                   *
 *  inputs:                                                                   *
 *      FuncName:                                                             *
 *          The name of the function at the python level.                     *
 *      CName:                                                                *
 *          The name of the function at the C level with the data type        *
 *          excluded. For example, if you wanted to wrap the Sinc_Double      *
 *          function defined in special_functions/, CName whould be Sinc.     *
 *  NOTE:                                                                     *
 *      Functions built using this macro MUST require three variables of the  *
 *      following types:                                                      *
 *          rho:                                                              *
 *              A numpy array or list of real numbers, a float, or an int.    *
 *          a:                                                                *
 *              A positive real number.                                       *
 *          F:                                                                *
 *              A positive real number.                                       *
 *      Ideally, only use this for the diffraction modeling functions which   *
 *      require three inputs. rho should be ring radius, a should be a        *
 *      specific radii, and F is the fresnel scale. rho, a, and F should all  *
 *      have the same units, preferably kilometers (km).                      */
#define RSS_RINGOCCSDiffModelingThreeVars(FuncName, CName)                     \
static PyObject *FuncName(PyObject *self, PyObject *args)                      \
{                                                                              \
    /*  We'll need output and capsule for safely creating the output array   */\
    /*  and ensuring we don't have a memory leak. rho is the input variable. */\
    PyObject *output, *capsule,  *rho;                                         \
                                                                               \
    /*  a and F are the input real numbers. Currently the code only          */\
    /*  accepts floats and ints, and does not deal with fancy Python objects.*/\
    double a, F;                                                               \
                                                                               \
    /*  Variable for the size of the input array which we'll retrieve later. */\
    long dim;                                                                  \
                                                                               \
    /*  Try to extract the values passed from the user, return error if fail.*/\
    if (!PyArg_ParseTuple(args, "Odd", &rho, &a, &F))                          \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rCould not parse inputs. Legal inputs are:\n"                    \
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"       \
            "\r\ta:     Positive constant (Float)\n"                           \
            "\r\tF      Positive constant (Float)\n\n"                         \
            "\rNotes:\n"                                                       \
            "\r\trho must be a non-empty one dimensional numpy array,\n"       \
            "\r\ta float, or an int.\n", VarToString(FuncName)                 \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  We're modeling in polar coordinates, so radius values should not be  */\
    /*  negative. Check this and raise an error if a < 0.                    */\
    if (a < 0.0)                                                               \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInner radius is not positive. (i.e. a<0)\n",                    \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  The Fresnel scale should never be zero or negative. Check.           */\
    else if (F <= 0.0)                                                         \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rFresnel scale is not positive (i.e. F<=0).\n",                  \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  If the user supplied a float or int, simply pass the inputs to the   */\
    /*  C function and return the value back to Python.                      */\
    if (PyLong_Check(rho) || PyFloat_Check(rho))                               \
    {                                                                          \
        /*  Convert a python float or int into a C double.                   */\
        double val = PyFloat_AsDouble(rho);                                    \
                                                                               \
        /*  The function returns a complex variable, so declare one.         */\
        complex double out;                                                    \
                                                                               \
        /*  Now get the computed value and convert the C complex double to a */\
        /*  Python complex value and return. The ## symbol tells the         */\
        /*  preprocessor to concatenate the words for the compiler. So if we */\
        /*  pass Gap_Diffraction as CName, the compiler will see             */\
        /*      out = Gap_Diffraction_Double(val, a, b, F);                  */\
        /*  which is exactly what we want.                                   */\
        out = CName##_Double(val, a, F);                                       \
                                                                               \
        /*  To build a Python complex value we need to pass the real and     */\
        /*  imaginary parts of out into PyComplex_FromDoubles. We can easily */\
        /*  access these values using creal and cimag, which are functions   */\
        /*  found in complex.h.                                              */\
        return PyComplex_FromDoubles(creal(out), cimag(out));                  \
    }                                                                          \
                                                                               \
    /*  As of v1.3 we allow lists to be passed to these routines.            */\
    else if (PyList_Check(rho))                                                \
    {                                                                          \
        /*  If the user passed a list, we'll need to return one. Create a    */\
        /*  variable for indexing over the input list.                       */\
        long i;                                                                \
                                                                               \
        /*  Create another variable for indexing over the list. This will be */\
        /*  the object corresponding to the value of the index i.            */\
        PyObject *ith_item;                                                    \
                                                                               \
        /*  And lastly, a double for storing the input value and a complex   */\
        /*  double for storing the output value.                             */\
        double val;                                                            \
        double complex out;                                                    \
                                                                               \
        /*  Get the number of elements in the list.                          */\
        dim = PyList_Size(rho);                                                \
                                                                               \
        /*  If the input list is empty, return with an error.                */\
        if (dim == 0)                                                          \
        {                                                                      \
            PyErr_Format(                                                      \
                PyExc_TypeError,                                               \
                "\n\rError Encountered: rss_ringoccs\n"                        \
                "\r\tdiffrec.ringlet_diffraction\n\n"                          \
                "\rInput list is empty.\n", VarToString(FuncName)              \
            );                                                                 \
            return NULL;                                                       \
        }                                                                      \
                                                                               \
        /*  Create a new list to be returned to the user.                    */\
        output = PyList_New(dim);                                              \
                                                                               \
        /*  Loop over the entries of the list, examine values for error, and */\
        /*  try to perform computations if legal inputs are given.           */\
        for (i=0; i<dim; ++i)                                                  \
        {                                                                      \
            /*  Extract the python object in the ith slot of the list.       */\
            ith_item = PyList_GET_ITEM(rho, i);                                \
                                                                               \
            /*  The list should be homogeneous with real numbers only. Check */\
            /*  each item and raise an error otherwise.                      */\
            if (!PyFloat_Check(ith_item) && !PyLong_Check(ith_item))           \
            {                                                                  \
                PyErr_Format(PyExc_TypeError,                                  \
                             "\n\rError Encountered: rss_ringoccs\n"           \
                             "\r\tdiffrec.%s\n\n"                              \
                             "\rInput list must contain real numbers only.\n", \
                             VarToString(FuncName));                           \
                return NULL;                                                   \
            }                                                                  \
                                                                               \
            /*  Convert the python float to a C double.                      */\
            val = PyFloat_AsDouble(ith_item);                                  \
                                                                               \
            /*  Compute the value, set the item in the list, and move on to  */\
            /*  to the next one. As explained before, ## concatenates words  */\
            /*  so CName##_Double compiles as CName_Double.                  */\
            out = CName##_Double(val, a, F);                                   \
            PyList_SET_ITEM(output, i,                                         \
                            PyComplex_FromDoubles(creal(out), cimag(out)));    \
        }                                                                      \
                                                                               \
        /*  The list created without errors, so just return it.              */\
        return output;                                                         \
    }                                                                          \
                                                                               \
    /*  If the input was not a list or a number, it must be a numpy array.   */\
    /*  Check this and return error otherwise.                               */\
    else if (!PyArray_Check(rho))                                              \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rCould not parse inputs. Legal inputs are:\n"                    \
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"       \
            "\r\ta:     Positive constant (Float)\n"                           \
            "\r\tb:     Positive constant (Float) greater than a\n"            \
            "\r\tF      Positive constant (Float)\n\n"                         \
            "\rNotes:\n"                                                       \
            "\r\trho must be a non-empty one dimensional numpy array,\n"       \
            "\r\ta float, or an int.\n", VarToString(FuncName)                 \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  Check the input array to make sure it is valid. It is important to   */\
    /*  check the number of dimensions BEFORE setting the dim variable. If   */\
    /*  the input array is empty, PyArray_DIMS(rho) will be a NULL pointer.  */\
    /*  Trying to access it with dim = PyArray_DIMS(rho)[0] will create a    */\
    /*  segmentation fault, crashing the Python interpreter.                 */\
    if (PyArray_NDIM((PyArrayObject *)rho) != 1)                               \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInput numpy array is not one-dimensional.\n",                   \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  Useful information about the data. The numpy API functions want a    */\
    /*  PyArrayObject pointer, so we need to cast rho as this type.          */\
    int typenum = PyArray_TYPE((PyArrayObject *)rho);                          \
    void *data  = PyArray_DATA((PyArrayObject *)rho);                          \
    dim         = PyArray_DIMS((PyArrayObject *)rho)[0];                       \
                                                                               \
    /*  This variable is for the data type of the output variable. It should */\
    /*  the complex valued equivalent of the input. So NPY_FLOAT should      */\
    /*  return NPY_CLOAT. All integer types will return NPY_CDOUBLE.         */\
    /*  If the input array was an int-like array (long, short, char, etc)    */\
    /*  then the output will be complex double by default/                   */\
    int outtype;                                                               \
                                                                               \
    /*  This void pointer will point to the output data we'll create later.  */\
    void *y;                                                                   \
                                                                               \
    /*  If the input array is empty, return with an error.                   */\
    if (dim == 0)                                                              \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInput numpy array is empty.\n", VarToString(FuncName)           \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  For float, double, and long double precision (which corresponds to   */\
    /*  numpy.float32, numpy.float64, and numpy.float128, respectively) we   */\
    /*  have functions which compute at these respective precisions. Pass    */\
    /*  the values into the _array_from_three macro which will loop over the */\
    /*  data and create a new output pointer, which y will then point to.    */\
    if (typenum == NPY_FLOAT)                                                  \
    {                                                                          \
        _array_from_three(data, a, F, y, dim, float,                           \
                          complex float, CName##_Float);                       \
        outtype = NPY_CFLOAT;                                                  \
    }                                                                          \
    else if (typenum == NPY_DOUBLE)                                            \
    {                                                                          \
        _array_from_three(data, a, F, y, dim, double,                          \
                          complex double, CName##_Double);                     \
        outtype = NPY_CDOUBLE;                                                 \
    }                                                                          \
    else if (typenum == NPY_LONGDOUBLE)                                        \
    {                                                                          \
        _array_from_three(data, a, F, y, dim, long double,                     \
                         long double, CName##_Long_Double);                    \
        outtype = NPY_CLONGDOUBLE;                                             \
    }                                                                          \
                                                                               \
    /*  For all other types we try to convert to double and compute.         */\
    else                                                                       \
    {                                                                          \
                                                                               \
        /*  Try to convert the input numpy array to double and compute.      */\
        PyObject *newrho = PyArray_FromObject(rho, NPY_DOUBLE, 1, 1);          \
                                                                               \
        /*  If PyArray_FromObject failed, newrho should be NULL. Check this. */\
        if (!(newrho))                                                         \
        {                                                                      \
            PyErr_Format(PyExc_TypeError,                                      \
                        "\n\rError Encountered: rss_ringoccs\n"                \
                        "\r\tdiffrec.%s\n\n"                                   \
                        "\rInvalid data type for input array. Input should"    \
                        "\n\rbe a 1-dimensional array of real numbers.\n",     \
                        VarToString(FuncName));                                \
            return NULL;                                                       \
        }                                                                      \
        else                                                                   \
            _array_from_three(data, a, F, y, dim, double,                      \
                              complex double, CName##_Double);                 \
                                                                               \
        /*  Retrieve a pointer to the data inside of newrho and set outtype. */\
        data = PyArray_DATA((PyArrayObject *)newrho);                          \
        outtype = NPY_CDOUBLE;                                                 \
    }                                                                          \
                                                                               \
    /*  Create a numpy array from the newly computed data and set a capsule  */\
    /*  for it. This ensures when the corresponding Python variable is       */\
    /*  deleted or removed the memory allocated to the respective C pointer  */\
    /*  is freed. Skipping this will result in memory leaks!                 */\
    output  = PyArray_SimpleNewFromData(1, &dim, outtype, y);                  \
    capsule = PyCapsule_New(y, NULL, capsule_cleanup);                         \
                                                                               \
    /*  This frees the variable at the Python level once it's destroyed.     */\
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);                   \
                                                                               \
    /*  Return the results to Python.                                     */   \
    return Py_BuildValue("N", output);                                         \
}

/*  Same as RSS_RINGOCCSDiffModelingFourVars but when the output is real and  *
 *  not complex. It has the following inputs:                                 *
 *      FuncName:                                                             *
 *          The name of the function at the python level.                     *
 *      CName:                                                                *
 *          The name of the function at the C level with the data type        *
 *          excluded. For example, if you wanted to wrap the Sinc_Double      *
 *          function defined in special_functions/, CName whould be Sinc.     *
 *  NOTE:                                                                     *
 *      Functions built using this macro MUST require four variables of the   *
 *      following types:                                                      *
 *          rho:                                                              *
 *              A numpy array or list of real numbers, a float, or an int.    *
 *          a:                                                                *
 *              A positive real number.                                       *
 *          b:                                                                *
 *              A positive real number GREATER than a.                        *
 *          F:                                                                *
 *              A positive real number.                                       *
 *      Ideally, only use this for the diffraction modeling functions which   *
 *      require four inputs and return reals, such as the phase computing     *
 *      routines. rho should be ring radius, a and b should be specific radii,*
 *      and F is the fresnel scale. rho, a, b, and F should all have the same *
 *      units, preferably kilometers (km). Output is usually radians.         */
#define RSS_RINGOCCSDiffModelingFourVarsREAL(FuncName, CName)                  \
static PyObject *FuncName(PyObject *self, PyObject *args)                      \
{                                                                              \
    /*  We'll need output and capsule for safely creating the output array   */\
    /*  and ensuring we don't have a memory leak. rho is the input variable. */\
    PyObject *output, *capsule,  *rho;                                         \
                                                                               \
    /*  a, b, and F are the input real numbers. Currently the code only      */\
    /*  accepts floats and ints, and does not deal with fancy Python objects.*/\
    double a, b, F;                                                            \
                                                                               \
    /*  Variable for the size of the input array which we'll retrieve later. */\
    long dim;                                                                  \
                                                                               \
    /*  Try to extract the values passed from the user, return error if fail.*/\
    if (!PyArg_ParseTuple(args, "Oddd", &rho, &a, &b, &F))                     \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rCould not parse inputs. Legal inputs are:\n"                    \
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"       \
            "\r\ta:     Positive constant (Float)\n"                           \
            "\r\tb:     Positive constant (Float) greater than a\n"            \
            "\r\tF      Positive constant (Float)\n\n"                         \
            "\rNotes:\n"                                                       \
            "\r\trho must be a non-empty one dimensional numpy array,\n"       \
            "\r\ta float, or an int.\n", VarToString(FuncName)                 \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  For all functions using this macro, a <= b. Check this.              */\
    if (a >= b)                                                                \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_ValueError,                                                  \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInner radius is not less than outer radius (i.e. a >= b).\n",   \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  We're modeling in polar coordinates, so radius values should not be  */\
    /*  negative. Check this and raise an error if a < 0.                    */\
    else if (a < 0.0)                                                          \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInner radius is not positive. (i.e. a<0)\n",                    \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  The Fresnel scale should never be zero or negative. Check.           */\
    else if (F <= 0.0)                                                         \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rFresnel scale is not positive (i.e. F<=0).\n",                  \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  If the user supplied a float or int, simply pass the inputs to the   */\
    /*  C function and return the value back to Python.                      */\
    if (PyLong_Check(rho) || PyFloat_Check(rho))                               \
    {                                                                          \
        /*  Convert a python float or int into a C double.                   */\
        double val = PyFloat_AsDouble(rho);                                    \
                                                                               \
        /*  Now get the computed value and convert the C complex double to a */\
        /*  Python complex value and return. The ## symbol tells the         */\
        /*  preprocessor to concatenate the words for the compiler. So if we */\
        /*  pass Gap_Diffraction as CName, the compiler will see             */\
        /*      out = Gap_Diffraction_Double(val, a, b, F);                  */\
        /*  which is exactly what we want.                                   */\
        val = CName##_Double(val, a, b, F);                                    \
                                                                               \
        /*  Now convert the C double back into a Python float and return.    */\
        return PyFloat_FromDouble(val);                                        \
    }                                                                          \
                                                                               \
    /*  As of v1.3 we allow lists to be passed to these routines.            */\
    else if (PyList_Check(rho))                                                \
    {                                                                          \
        /*  If the user passed a list, we'll need to return one. Create a    */\
        /*  variable for indexing over the input list.                       */\
        long i;                                                                \
                                                                               \
        /*  Create another variable for indexing over the list. This will be */\
        /*  the object corresponding to the value of the index i.            */\
        PyObject *ith_item;                                                    \
                                                                               \
        /*  And lastly, a double for storing the input value.                */\
        double val;                                                            \
                                                                               \
        /*  Get the number of elements in the list.                          */\
        dim = PyList_Size(rho);                                                \
                                                                               \
        /*  If the input list is empty, return with an error.                */\
        if (dim == 0)                                                          \
        {                                                                      \
            PyErr_Format(                                                      \
                PyExc_TypeError,                                               \
                "\n\rError Encountered: rss_ringoccs\n"                        \
                "\r\tdiffrec.%s\n\n"                                           \
                "\rInput list is empty.\n", VarToString(FuncName)              \
            );                                                                 \
            return NULL;                                                       \
        }                                                                      \
                                                                               \
        /*  Create a new list to be returned to the user.                    */\
        output = PyList_New(dim);                                              \
                                                                               \
        /*  Loop over the entries of the list, examine values for error, and */\
        /*  try to perform computations if legal inputs are given.           */\
        for (i=0; i<dim; ++i)                                                  \
        {                                                                      \
            /*  Extract the python object in the ith slot of the list.       */\
            ith_item = PyList_GET_ITEM(rho, i);                                \
                                                                               \
            /*  The list should be homogeneous with real numbers only. Check */\
            /*  each item and raise an error otherwise.                      */\
            if (!PyFloat_Check(ith_item) && !PyLong_Check(ith_item))           \
            {                                                                  \
                PyErr_Format(PyExc_TypeError,                                  \
                             "\n\rError Encountered: rss_ringoccs\n"           \
                             "\r\tdiffrec.%s\n\n"                              \
                             "\rInput list must contain real numbers only.\n", \
                             VarToString(FuncName));                           \
                return NULL;                                                   \
            }                                                                  \
                                                                               \
            /*  Convert the python float to a C double.                      */\
            val = PyFloat_AsDouble(ith_item);                                  \
                                                                               \
            /*  Compute the value, set the item in the list, and move on to  */\
            /*  to the next one. As explained before, ## concatenates words  */\
            /*  so CName##_Double compiles as CName_Double.                  */\
            val = CName##_Double(val, a, b, F);                                \
            PyList_SET_ITEM(output, i, PyFloat_FromDouble(val));               \
        }                                                                      \
                                                                               \
        /*  The list created without errors, so just return it.              */\
        return output;                                                         \
    }                                                                          \
                                                                               \
    /*  If the input was not a list or a number, it must be a numpy array.   */\
    /*  Check this and return error otherwise.                               */\
    else if (!PyArray_Check(rho))                                              \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rCould not parse inputs. Legal inputs are:\n"                    \
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"       \
            "\r\ta:     Positive constant (Float)\n"                           \
            "\r\tb:     Positive constant (Float) greater than a\n"            \
            "\r\tF      Positive constant (Float)\n\n"                         \
            "\rNotes:\n"                                                       \
            "\r\trho must be a non-empty one dimensional numpy array,\n"       \
            "\r\ta float, or an int.\n", VarToString(FuncName)                 \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  Check the input array to make sure it is valid. It is important to   */\
    /*  check the number of dimensions BEFORE setting the dim variable. If   */\
    /*  the input array is empty, PyArray_DIMS(rho) will be a NULL pointer.  */\
    /*  Trying to access it with dim = PyArray_DIMS(rho)[0] will create a    */\
    /*  segmentation fault, crashing the Python interpreter.                 */\
    if (PyArray_NDIM((PyArrayObject *)rho) != 1)                               \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInput numpy array is not one-dimensional.\n",                   \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  Useful information about the data. The numpy API functions want a    */\
    /*  PyArrayObject pointer, so we need to cast rho as this type.          */\
    int typenum = PyArray_TYPE((PyArrayObject *)rho);                          \
    void *data  = PyArray_DATA((PyArrayObject *)rho);                          \
    dim         = PyArray_DIMS((PyArrayObject *)rho)[0];                       \
                                                                               \
    /*  This void pointer will point to the output data we'll create later.  */\
    void *y;                                                                   \
                                                                               \
    /*  If the input array is empty, return with an error.                   */\
    if (dim == 0)                                                              \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInput numpy array is empty.\n", VarToString(FuncName)           \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  For float, double, and long double precision (which corresponds to   */\
    /*  numpy.float32, numpy.float64, and numpy.float128, respectively) we   */\
    /*  have functions which compute at these respective precisions. Pass    */\
    /*  the values into the _array_from_four macro which will loop over the  */\
    /*  data and create a new output pointer, which y will then point to.    */\
    if (typenum == NPY_FLOAT)                                                  \
        _array_from_four(data, a, b, F, y, dim, float, float, CName##_Float);  \
    else if (typenum == NPY_DOUBLE)                                            \
        _array_from_four(data, a, b, F, y, dim, double,                        \
                         double, CName##_Double);                              \
    else if (typenum == NPY_LONGDOUBLE)                                        \
        _array_from_four(data, a, b, F, y, dim, long double,                   \
                         long double, CName##_Long_Double);                    \
                                                                               \
    /*  For all other types we try to convert to double and compute.         */\
    else                                                                       \
    {                                                                          \
                                                                               \
        /*  Try to convert the input numpy array to double and compute.      */\
        PyObject *newrho = PyArray_FromObject(rho, NPY_DOUBLE, 1, 1);          \
                                                                               \
        /*  If PyArray_FromObject failed, newrho should be NULL. Check this. */\
        if (!(newrho))                                                         \
        {                                                                      \
            PyErr_Format(PyExc_TypeError,                                      \
                        "\n\rError Encountered: rss_ringoccs\n"                \
                        "\r\tdiffrec.%s\n\n"                                   \
                        "\rInvalid data type for input array. Input should"    \
                        "\n\rbe a 1-dimensional array of real numbers.\n",     \
                        VarToString(FuncName));                                \
            return NULL;                                                       \
        }                                                                      \
        else                                                                   \
            _array_from_four(data, a, b, F, y, dim, double,                    \
                             double, CName##_Double);                          \
                                                                               \
        /*  Retrieve a pointer to the data inside of newrho.                 */\
        data = PyArray_DATA((PyArrayObject *)newrho);                          \
    }                                                                          \
                                                                               \
    /*  Create a numpy array from the newly computed data and set a capsule  */\
    /*  for it. This ensures when the corresponding Python variable is       */\
    /*  deleted or removed the memory allocated to the respective C pointer  */\
    /*  is freed. Skipping this will result in memory leaks!                 */\
    output  = PyArray_SimpleNewFromData(1, &dim, typenum, y);                  \
    capsule = PyCapsule_New(y, NULL, capsule_cleanup);                         \
                                                                               \
    /*  This frees the variable at the Python level once it's destroyed.     */\
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);                   \
                                                                               \
    /*  Return the results to Python.                                     */   \
    return Py_BuildValue("N", output);                                         \
}

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
 *      if not hasattr(dlpin, attr_name):                                     *
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
#define RSS_RINGOCCSCheckArray(ptr, dlpin, attr_name, tmp, arr, length)        \
if (!PyObject_HasAttrString(dlpin, VarToString(attr_name)))                    \
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
    tmp = PyObject_GetAttrString(dlpin, VarToString(attr_name));               \
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
else if (PyArray_DIMS((PyArrayObject *)arr)[0] != length)                      \
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
#define RSS_RINGOCCSSetArray(self, dlp, VarName, tmp, typenum, start, length)  \
PyObject *PyArray_##VarName;                                                   \
PyObject *PyCapsule_##VarName;                                                 \
double *VarName = (double *)malloc(sizeof(double) * length);                   \
                                                                               \
for (i=0; i<length; ++i)                                                       \
    VarName[i] = dlp.VarName[i+start];                                         \
                                                                               \
PyArray_##VarName = PyArray_SimpleNewFromData(1, &length, typenum,             \
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
if (Min_Double(varname, arrsize) <= 0.0)                                       \
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
if ((Min_Double(varname, arrsize) < -TWO_PI) ||                                \
    (Max_Double(varname, arrsize) > TWO_PI))                                   \
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

#define RSS_RINGOCCSPsitypeError                                               \
    PyErr_Format(                                                              \
        PyExc_ValueError,                                                      \
        "\n\rError Encountered: rss_ringoccs\n"                                \
        "\r\tdiffrec.DiffractionCorrection\n\n"                                \
        "\r\tIllegal string for psitype. Allowed strings are:\n"               \
        "\r\t\tnewton:     Newton-Raphson method\n"                            \
        "\r\t\tfresnel:    Quadratic Fresnel approximation\n"                  \
        "\r\t\tfresneln:   Legendre polynomial approximation with 1<n<256\n"   \
    );

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

/*  This function acts like a Python dictionary and returns the correct       *
 *  normalized equivalent width for the corresponding input string. If an     *
 *  invalid string is given this returns -1. The user of this function then   *
 *  has the responsibility of handling this value appropriately, no errors    *
 *  are raised directly from this routine.                                    *
 *                                                                            *
 *  A small note on strcmp (string-compare). This does not simply return 0 or *
 *  1 depending on whether the two inputs are identical, but rather can       *
 *  return negative, positive, or zero, depending on whether one string is    *
 *  "greater" than the other. For identical strings this returns 0.           *
 *  Hence DO NOT write something like:                                        *
 *      if (strcmp(string, same_string))                                      *
 *          do stuff;                                                         *
 *  since strcmp(string, same_string) will return 0, meaning the if-statement *
 *  is false. Instead, do the following:                                      *
 *      if (strcmp(string, same_string) == 0)                                 *
 *          do stuff;                                                         */
static void GetNormeqFromString(char *wtype, double *norm_eq)
{
    if (strcmp(wtype, "rect") == 0)
        *norm_eq = RectNormEQ;
    else if (strcmp(wtype, "coss") == 0)
        *norm_eq = CossNormEQ;
    else if (strcmp(wtype, "kb20") == 0)
        *norm_eq = KB20NormEQ;
    else if (strcmp(wtype, "kb25") == 0)
        *norm_eq = KB25NormEQ;
    else if (strcmp(wtype, "kb35") == 0)
        *norm_eq = KB35NormEQ;
    else if (strcmp(wtype, "kbmd20") == 0)
        *norm_eq = KBMD20NormEQ;
    else if (strcmp(wtype, "kbmd25") == 0)
        *norm_eq = KBMD25NormEQ;
    else if (strcmp(wtype, "kbmd35") == 0)
        *norm_eq = KBMD35NormEQ;
    else
        *norm_eq = -1.0;
}

/*  This function provides a plethora of strings for the user to provide for  *
 *  the sake of convenience in specifying features of the rings of Saturn.    *
 *  This takes in a string and a pointer to the range variable, and on        *
 *  success it changes the values rng_list points to and then returns 1. On   *
 *  failure it return 0. Hence you may use this function like:                *
 *      if (!GetRangeFromString(my_string, my_range))                         *
 *          RaiseAnError;                                                     *
 *  All of the features correspond to Saturn. The user can append their own   *
 *  features for other planets, if desired, but MAKE SURE YOU ALTER THE ERROR *
 *  CHECKING ROUTINE found in the Diffrec_init function. If you do not, even  *
 *  if you've added a new string to this function, Diffrec_init will refuse   *
 *  to accept your new string and raise an error.                             */
static int GetRangeFromString(char *range, double *rng_list)
{
    if (strcmp(range, "all") == 0)
    {
        rng_list[0] = 1.0;
        rng_list[1] = 400000.0;
    }
    else if ((strcmp(range, "besselbarnard") == 0) ||
             (strcmp(range, "bessel-barnard") == 0))
    {
        rng_list[0] = 120210.0;
        rng_list[1] = 120330.0;
    }

    /*  This region is only interesting for low-inclination occultations,     *
     *  such as (Cassini) Rev133 and revs of similar occultation number.      */
    else if (strcmp(range, "cringripples") == 0)
    {
        rng_list[0] = 77690.0;
        rng_list[1] = 77760.0;
    }

    /*  Excellent test-case for most occultations since the Encke gap acts    *
     *  like a near perfect square well. Examining this region is worthwhile. */
    else if ((strcmp(range, "encke") == 0) || (strcmp(range, "enckegap") == 0))
    {
        rng_list[0] = 132900.0;
        rng_list[1] = 134200.0;
    }
    else if ((strcmp(range, "herschel") == 0) ||
             (strcmp(range, "herschelgap") == 0))
    {
        rng_list[0] = 118100.0;
        rng_list[1] = 118380.0;
    }
    else if ((strcmp(range, "huygens") == 0) ||
             (strcmp(range, "huygensringlet") == 0))
    {
        rng_list[0] = 117650.0;
        rng_list[1] = 117950.0;
    }
    else if (strcmp(range, "janusepimetheus") == 0)
    {
        rng_list[0] = 96200.0;
        rng_list[1] = 96800.0;
    }
    else if ((strcmp(range, "jeffreys") == 0) ||
             (strcmp(range, "jeffreysgap") == 0))
    {
        rng_list[0] = 118900.0;
        rng_list[1] = 119000.0;
    }
    else if ((strcmp(range, "kuiper") == 0) ||
             (strcmp(range, "kuipergap") == 0))
    {
        rng_list[0] = 119300.0;
        rng_list[1] = 119500.0;
    }

    /*  The go-to test-case for Team-Cassini 2017-2019. An excellent gap to   *
     *  examine and test routines on. The excellent geometry of Rev007 means  *
     *  one can get very accurate results on the Maxwell ringlet using any of *
     *  the reconstruction methods.                                           */
    else if ((strcmp(range, "maxwell") == 0) ||
             (strcmp(range, "maxwellringlet") == 0))
    {
        rng_list[0] = 87410.0;
        rng_list[1] = 87610.0;
    }
    else if ((strcmp(range, "russell") == 0) ||
             (strcmp(range, "russellgap") == 0))
    {
        rng_list[0] = 118550.0;
        rng_list[1] = 118660.0;
    }

    /*  A very steep gap, excellent for testing the resolution of the         *
     *  reconstruction.                                                       */
    else if ((strcmp(range, "titan") == 0) ||
             (strcmp(range, "titanringlet") == 0))
    {
        rng_list[0] = 77870.0;
        rng_list[1] = 77930.0;
    }
    else return 0;

    return 1;
}

/*  This function frees the memory allocated to a pointer by malloc when the  *
 *  corresponding variable is destroyed at the Python level. Without this you *
 *  will have serious memory leaks, so do not remove!                         */
static void capsule_cleanup(PyObject *capsule)
{
    void *memory = PyCapsule_GetPointer(capsule, NULL);
    free(memory);
}

/*  Standard algorithm of time-complexity O(n) to take an array and reverse   *
 *  the order. This is equivalent to taking a numpy array arr and writing     *
 *  arr = arr[::-1] (in Python 3, at least). Since all of the main routines   *
 *  use double pointers, only a double version is provided.                   */
static void ReverseDoubleArray(double *arr, long arrsize)
{
    double val;
    long i;
    for(i=0; i<arrsize/2; i++)
    {
        val = arr[i];
        arr[i] = arr[arrsize-i-1];
        arr[arrsize-i-1] = val;
    }
}

/*  Deallocating function for the DiffractionCorrection class.                */
static void Diffrec_dealloc(PyDiffrecObj *self)
{
    Py_XDECREF(self->B_rad_vals);
    Py_XDECREF(self->D_km_vals);
    Py_XDECREF(self->F_km_vals);
    Py_XDECREF(self->T_hat_fwd_vals);
    Py_XDECREF(self->T_hat_vals);
    Py_XDECREF(self->T_vals);
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
    long k;

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
    self->fwd = false;
    self->use_fft = false;
    self->verbose = false;

    /*  Using the bfac guarantees accurate window sizes in the case of a poor *
     *  Allen deviation. Window normalization is also recommended since the   *
     *  integral is scaled by the width of the window, and hence for small    *
     *  window sizes the result might return close to zero.                   */
    self->bfac = true;
    self->use_norm = true;

    /*  The default sigma value is the one for Cassini.                       */
    self->sigma = 2.0e-13;

    /*  If res_factor was not set, set to 0.75. This value was specified by   *
     *  Essam Marouf as necessary to ensure the reconstruction matches the    *
     *  PDS results. No justification is known to me.                         */
    self->res_factor = 0.75;

    /*  The write files option allows the user to export the end result to a  *
     *  CSV. Default is set to false.                                         */
    self->write_file = false;

    /*  The default geometry assumes the rings are circular, so we set both   *
     *  the eccentricity and the periapse to zero.                            */
    self->ecc = 0.0;
    self->peri = 0.0;

    /*  By default, polynomial interpolation is off.                          */
    self->interp = 0;

    /*  Default polynomial perturbation is off.                               */
    perturb = NULL;

    /*  The perturb variable is a list and should contain exactly 5 elements. */
    #define PERTURB_LEN 5

    /*  If the range input is a list, it should have exactly two elements.    */
    #define RANGE_LEN 2

    /*  Useful numbers for checking the input psitype and wtype.              */
    #define FRESNEL_STR_LENGTH 7
    #define BASETEN 10

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
    DLPObj dlp;

    /*  Check that interp is a legal value.                                   */
    if ((0 <= self->interp) && (self->interp <=4))
        dlp.interp = self->interp;
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

        /*  dlp.perturb still needs to be set, so loop over and set to zero.  */
        for (i=0; i<PERTURB_LEN; ++i)
            dlp.perturb[i] = 0.0;

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
         *  to doubles, and store them in the dlp.perturb variable.           */
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
                dlp.perturb[i] = (double)PyLong_AsLong(next);

            /*  Convert from Python float to C double with PyFloat_AsDouble.  */
            else if (PyFloat_Check(next))
                dlp.perturb[i] = PyFloat_AsDouble(next);

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

    /*  Run error checks on the input strings, including conversion to lower  *
     *  case. Note that this does NOT remove white space from the strings,    *
     *  additional white space will result in an invalid string. It is a bit  *
     *  of a pain to resize the strings and remove white space, so instead    *
     *  this program only accepts strings without additional spaces.          */

    /*  Get the size of the wtype and psitype strings. strlen is a C standard *
     *  library function that can be found in string.h.                       */
    long wtype_length   = strlen(self->wtype);
    long psitype_length = strlen(self->psitype);

    /*  Create variables to convert wtype and psitype to lowercase. Malloc is *
     *  a standard library function for allocating memory found in stdlib.h.  *
     *  we allocate +1 more than needed since strings in C are NULL           *
     *  terminated, and hence we need extra room for the NULL.                */
    dlp.wtype   = (char *)malloc(sizeof(char) * wtype_length + 1);
    dlp.psitype = (char *)malloc(sizeof(char) * wtype_length + 1);

    /*  Make a copy of the strings from python. Do NOT alter actual strings.  *
     *  strcpy is a C standard library function found in string.h.            */
    strcpy(dlp.wtype, self->wtype);
    strcpy(dlp.psitype, self->psitype);

    /*  Loop over the strings and convert to lower case. tolower is a C       *
     *  standard library function found in ctype.h, and NOT string.h. The     *
     *  reasons are unbeknownst to be, but knownst to others.                 */
    for (i=0; i<wtype_length; ++i)
        dlp.wtype[i] = tolower(self->wtype[i]);

    for (i=0; i<psitype_length; ++i)
        dlp.psitype[i] = tolower(self->psitype[i]);

    /*  Check that wtype is a legal value. Return ValueError otherwise.       *
     *  strcmp is found in string.h and can return positive, negative, or     *
     *  zero depending on if one string is "larger" than the other. It        *
     *  returns zero if the strings are identical, so DO NOT set:             *
     *      if (strcmp(string1, string2))                                     *
     *          do stuff;                                                     *
     *  since strcmp will return 0 if string1 = string2.                      */
    if (!(strcmp(dlp.wtype, "rect") == 0) &&
        !(strcmp(dlp.wtype, "coss") == 0) &&
        !(strcmp(dlp.wtype, "kb20") == 0) &&
        !(strcmp(dlp.wtype, "kb25") == 0) &&
        !(strcmp(dlp.wtype, "kb35") == 0) &&
        !(strcmp(dlp.wtype, "kbmd20") == 0) &&
        !(strcmp(dlp.wtype, "kbmd25") == 0) &&
        !(strcmp(dlp.wtype, "kbmd35") == 0))
    {
        PyErr_Format(
            PyExc_ValueError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r\tIllegal string for wtype. Allowed strings are:\n"
            "\r\t\trect:    Rectangular Window\n"
            "\r\t\tcoss:    Squared Cosine Window\n"
            "\r\t\tkb20:    Kaiser-Bessel with alpha=2.0 pi\n"
            "\r\t\tkb25:    Kaiser-Bessel with alpha=2.5 pi\n"
            "\r\t\tkb35:    Kaiser-Bessel with alpha=3.5 pi\n"
            "\r\t\tkbmd20:  Modified Kaiser-Bessel with alpha=2.0 pi\n"
            "\r\t\tkbmd25:  Modified Kaiser-Bessel with alpha=2.5 pi\n"
            "\r\t\tkbmd35:  Modified Kaiser-Bessel with alpha=3.5 pi\n"
        );
        return -1;
    }

    /*  Check that psitype is legal. The strstr is a C standard library       *
     *  function (from string.h) that checks if the second string is          *
     *  contained in the first. So strstr(dlp.psitype, "fresnel") is          *
     *  equivalent to the Python code ("fresnel" in dlp.psitype)              */
    if (!(strcmp(dlp.psitype, "newton") == 0) &&
         (strstr(dlp.psitype, "fresnel") == NULL))
    {
        /*  The macro RSS_RINGOCCSPsitypeError already has a semi-colon at    *
         *  the end so we don't need to add a new one here.                   */
        RSS_RINGOCCSPsitypeError
        return -1;
    }

    /*  If psitype is just "fresnel", set order to 1. The code will           *
     *  automatically pass the dlp instance to the correct function.          */
    if (strcmp(dlp.psitype, "fresnel") == 0)
        dlp.order = 1;

    /*  If psitype is "newton", set order to 1. This will mean the dlp        *
     *  instance is passed to the Newton-Raphson function.                    */
    else if (strcmp(dlp.psitype, "newton") == 0)
        dlp.order = 0;

    /*  strncmp is a C standard library function that compares the first n    *
     *  elements of two strings. If the first seven elements of dlp.psitype   *
     *  are "fresnel", but the string is not exactly "fresnel", try to parse  *
     *  the rest of it and extract a value. For example, if                   *
     *  dlp.psitype = "fresnel4", try to extract the "4".                     */
    else if (strncmp(dlp.psitype, "fresnel", FRESNEL_STR_LENGTH) == 0)
    {
        /*  Create a variable equal to the last elements of dlp.psitype.      */
        const char *psitype_num = &dlp.psitype[FRESNEL_STR_LENGTH];

        /*  Try to convert the last elements to an integer using the strtol   *
         *  function, which is a part of the C standard library. It stands    *
         *  for "string to long". Equivalent of (int)"number" in Python.      */
        dlp.order = strtol(psitype_num, NULL, BASETEN);

        /*  strtol returns 0 if it couldn't convert string to number. If this *
         *  happens the user likely provided a faulty string. Raise an error. */
        if (dlp.order == 0)
        {
            RSS_RINGOCCSPsitypeError
            return -1;
        }
    }

    /*  Illegal psitype was passed, return with error.                        */
    else
    {
        RSS_RINGOCCSPsitypeError
        return -1;
    }

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
        dlp.ecc = self->ecc;

    /*  Check that peri lies in the range -2pi<peri<2pi.                      */
    if (fabs(self->peri)>TWO_PI)
    {
        RSS_RINGOCCSTwoPiValueError(DiffractionCorrection, peri);
        return -1;
    }
    else
        dlp.peri = self->peri;

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
        long rng_length = strlen(rng_string);

        /*  Convert the string to lower-case.                                 */
        for (i=0; i<rng_length; ++i)
            rng_string[i] = tolower(rng_string[i]);

        /*  Pass the string into the get range function. This function stores *
         *  the corresponding values in self->range and returns 1 on success, *
         *  and returns 0 on failure. If fails, return with error.            */
        if (!GetRangeFromString(rng_string, self->range))
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
        long rng_length = strlen(rng_string_cpy);

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
        if (!GetRangeFromString(rng_string, self->range))
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

    /*  If every passed, set dlp.rho_km_vals to point to the data inside arr. */
    dlp.rho_km_vals = (double *)PyArray_DATA((PyArrayObject *)arr);
    dlp.arr_size = PyArray_DIMS((PyArrayObject *)arr)[0];

    /*  Check that there's actually data to process.                          */
    if (dlp.arr_size<2)
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
    RSS_RINGOCCSCheckArray(dlp.rho_km_vals, DLPInst, rho_km_vals,
                           tmp, arr, dlp.arr_size);
    RSS_RINGOCCSCheckArray(dlp.B_rad_vals, DLPInst, B_rad_vals,
                           tmp, arr, dlp.arr_size);
    RSS_RINGOCCSCheckArray(dlp.phi_rad_vals, DLPInst, phi_rad_vals,
                           tmp, arr, dlp.arr_size);
    RSS_RINGOCCSCheckArray(dlp.D_km_vals, DLPInst, D_km_vals,
                           tmp, arr, dlp.arr_size);
    RSS_RINGOCCSCheckArray(dlp.f_sky_hz_vals, DLPInst, f_sky_hz_vals,
                           tmp, arr, dlp.arr_size);
    RSS_RINGOCCSCheckArray(dlp.rho_dot_kms_vals, DLPInst,
                           rho_dot_kms_vals, tmp, arr, dlp.arr_size);
    RSS_RINGOCCSCheckArray(dlp.t_oet_spm_vals, DLPInst,
                           t_oet_spm_vals, tmp, arr, dlp.arr_size);
    RSS_RINGOCCSCheckArray(dlp.t_ret_spm_vals, DLPInst,
                           t_ret_spm_vals, tmp, arr, dlp.arr_size);
    RSS_RINGOCCSCheckArray(dlp.t_set_spm_vals, DLPInst,
                           t_set_spm_vals, tmp, arr, dlp.arr_size);
    RSS_RINGOCCSCheckArray(dlp.rho_corr_pole_km_vals, DLPInst,
                           rho_corr_pole_km_vals,tmp, arr, dlp.arr_size);
    RSS_RINGOCCSCheckArray(dlp.rho_corr_timing_km_vals, DLPInst,
                           rho_corr_timing_km_vals, tmp, arr, dlp.arr_size);
    RSS_RINGOCCSCheckArray(dlp.raw_tau_threshold_vals, DLPInst,
                           raw_tau_threshold_vals, tmp, arr, dlp.arr_size);
    RSS_RINGOCCSCheckArray(dlp.phi_rl_rad_vals, DLPInst,
                           phi_rl_rad_vals, tmp, arr, dlp.arr_size);
    RSS_RINGOCCSCheckArray(dlp.p_norm_vals, DLPInst, p_norm_vals,
                           tmp, arr, dlp.arr_size);
    RSS_RINGOCCSCheckArray(dlp.phase_rad_vals, DLPInst,
                           phase_rad_vals, tmp, arr, dlp.arr_size);

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

    /*  The input phase needs to be negated due to mathematical conventions.  */
    for (i=0; i<dlp.arr_size; ++i)
        dlp.phase_rad_vals[i] = -dlp.phase_rad_vals[i];

    /*  Store the dlp history inside of the DiffractionCorrection class. This *
     *  tmp, Py_INCREF, Py_XDECREF method is recommended in the Python C-API  *
     *  documentation as a means of safely storing the variable.              */
    tmp = self->dathist;
    Py_INCREF(history);
    self->dathist = history;
    Py_XDECREF(tmp);

    RSS_RINGOCCSCheckArrayPositive(dlp.D_km_vals, dlp.arr_size);
    RSS_RINGOCCSCheckArrayPositive(dlp.rho_km_vals, dlp.arr_size);
    RSS_RINGOCCSCheckArrayPositive(dlp.f_sky_hz_vals, dlp.arr_size);
    RSS_RINGOCCSCheckArrayTwoPi(dlp.B_rad_vals, dlp.arr_size);
    RSS_RINGOCCSCheckArrayTwoPi(dlp.phi_rad_vals, dlp.arr_size);

    /*  Make sure the power isn't negative.                                   */
    if (Min_Double(dlp.p_norm_vals, dlp.arr_size) < 0.0)
    {
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rdlp.p_norm_vals contains non-positive values.\n"
        );
        return -1;
    }

    /*  Run an error check as to whether or not this is an egress, ingress,   *
     *  or chord occultation, and if the data can be processed accordingly.   */
    double min_val = Min_Double(dlp.rho_dot_kms_vals, dlp.arr_size);
    double max_val = Max_Double(dlp.rho_dot_kms_vals, dlp.arr_size);
    double dx_km   = dlp.rho_km_vals[1] - dlp.rho_km_vals[0];

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
        for(i=0; i<dlp.arr_size; ++i)
            dlp.rho_dot_kms_vals[i] = fabs(dlp.rho_dot_kms_vals[i]);
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
        ReverseDoubleArray(dlp.rho_km_vals, dlp.arr_size);
        ReverseDoubleArray(dlp.phi_rad_vals, dlp.arr_size);
        ReverseDoubleArray(dlp.B_rad_vals, dlp.arr_size);
        ReverseDoubleArray(dlp.D_km_vals, dlp.arr_size);
        ReverseDoubleArray(dlp.phase_rad_vals, dlp.arr_size);
        ReverseDoubleArray(dlp.p_norm_vals, dlp.arr_size);
        ReverseDoubleArray(dlp.f_sky_hz_vals, dlp.arr_size);
        ReverseDoubleArray(dlp.rho_dot_kms_vals, dlp.arr_size);
        for(i=0; i<dlp.arr_size; ++i)
            dlp.rho_dot_kms_vals[i] = fabs(dlp.rho_dot_kms_vals[i]);
        dx_km *= -1.0;
    }

    if (self->verbose)
        puts("\tComputing Necessary Variables...");

    double *lambda_sky = (double *)malloc(sizeof(double) * dlp.arr_size);
    double *mu_vals    = (double *)malloc(sizeof(double) * dlp.arr_size);

    /*  Allocate memory for the diffracted data and the fresnel scale.        */
    dlp.T_in = (complex double *)malloc(sizeof(complex double) * dlp.arr_size);
    dlp.F_km_vals = (double *)malloc(sizeof(double) * dlp.arr_size);

    for (i=0; i<dlp.arr_size; ++i)
    {
        /*  Compute the complex amplitude, T_hat_vals.                        */
        dlp.T_in[i] = (cos(dlp.phase_rad_vals[i]) +
                       _Complex_I*sin(dlp.phase_rad_vals[i]))
                      *sqrt(dlp.p_norm_vals[i]);

        /*  Compute the wavelength lambda and scale-factor mu.                */
        lambda_sky[i] = Frequency_To_Wavelength_Double(dlp.f_sky_hz_vals[i]);
        mu_vals[i]    = sin(fabs(dlp.B_rad_vals[i]));

        /*  And finally, compute the Fresnel scale.                           */
        dlp.F_km_vals[i]  = Fresnel_Scale_Double(lambda_sky[i],
                                                 dlp.D_km_vals[i],
                                                 dlp.phi_rad_vals[i],
                                                 dlp.B_rad_vals[i]);
    }

    /*  Get the normalized equivalent width of the window type.               */
    GetNormeqFromString(dlp.wtype, &self->norm_eq);

    /*  Use calloc to both allocate memory for dlp.w_km_vals (like malloc)    *
     *  and initialize the data to zero (unlike malloc). This is similar to   *
     *  numpy.zeros(dlp.arr_size) in Python.                                  */
    dlp.w_km_vals = (double *)calloc(dlp.arr_size, sizeof(double));

    /*  Declare long pointer-to-pointer which stores the indices where        *
     *  F_km_vals is non-zero in the first slot (Prange[0]), and the size of  *
     *  this array in the second (*Prange[1]).                                */
    long **Prange;
    long *Prange_Index;
    long Prange_Size;
    double w_fac;

    if (self->bfac)
    {
        w_fac = self->norm_eq;
        double omega;
        double *alpha  = (double *)malloc(sizeof(double) * dlp.arr_size);
        double *P_vals = (double *)malloc(sizeof(double) * dlp.arr_size);

        for(i=0; i<dlp.arr_size; ++i)
        {
            omega      = TWO_PI * dlp.f_sky_hz_vals[i];
            alpha[i]   = omega * self->sigma;
            alpha[i]  *= alpha[i] * 0.5 / dlp.rho_dot_kms_vals[i];
            P_vals[i]  = self->res/(alpha[i]*dlp.F_km_vals[i]*dlp.F_km_vals[i]);
        }

        Prange = Where_Greater_Double(P_vals, dlp.arr_size, 1.0);
        Prange_Index = Prange[0];
        Prange_Size  = *Prange[1];

        for (i=0; i<Prange_Size; ++i)
            dlp.w_km_vals[Prange_Index[i]] = w_fac * 
                Resolution_Inverse_Double(P_vals[Prange_Index[i]]) / alpha[i];

        free(P_vals);
        free(alpha);
    }
    else
    {
        w_fac = self->norm_eq/self->res;

        for (i=0; i<dlp.arr_size; ++i)
            dlp.w_km_vals[i] = 2.0*dlp.F_km_vals[i]*dlp.F_km_vals[i]*w_fac;

        Prange = Where_Greater_Double(dlp.F_km_vals, dlp.arr_size, 0.0);
        Prange_Index = Prange[0];
        Prange_Size  = *Prange[1];
    }

    double *rho_legal = (double *)malloc(sizeof(double) * Prange_Size);

    for(i=0; i<Prange_Size; ++i)
        rho_legal[i] = dlp.rho_km_vals[Prange_Index[i]] -
                       0.5*dlp.w_km_vals[Prange_Index[i]];

    max_val = Max_Double(rho_legal, Prange_Size);

    for(i=0; i<Prange_Size; ++i)
        rho_legal[i] = dlp.rho_km_vals[Prange_Index[i]] +
                       0.5*dlp.w_km_vals[Prange_Index[i]];
    min_val = Min_Double(rho_legal, Prange_Size);

    long **wrange = Where_LesserGreater_Double(dlp.rho_km_vals, dlp.arr_size,
                                                min_val, max_val);
    long *wrange_Index = wrange[0];
    long wrange_Size = *wrange[1];

    if (wrange_Size == 0)
    {
        PyErr_Format(
            PyExc_ValueError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r\tThe window width is too large to reconstruct anything.\n"
            "\r\t\tMinimum Available Radius:         %f\n"
            "\r\t\tMaximum Available Radius:         %f\n"
            "\r\t\tMinimum Required Window Width:    %f\n"
            "\r\t\tMaximum Required Window Width:    %f\n",
            min_val, max_val, Min_Double(dlp.w_km_vals, dlp.arr_size),
            Max_Double(dlp.w_km_vals, dlp.arr_size)
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
        self->start = wrange_Index[0];
        self->finish = wrange_Index[wrange_Size-1];
        self->n_used = self->finish-self->start;
    }

    dlp.start  = self->start;
    dlp.n_used = self->n_used;

    /*  For the first inverse calculatiion, dlp.use_fwd must be set to false. */
    dlp.use_fwd   = false;

    if (self->use_norm)
        dlp.use_norm = true;

    if (self->use_fft)
        dlp.use_fft = true;

    dlp.kd_vals = (double *)malloc(sizeof(double) * dlp.arr_size);
    for(i=0; i<dlp.arr_size; ++i)
        dlp.kd_vals[i] = Wavelength_To_Wavenumber_Double(lambda_sky[i]) *
                     dlp.D_km_vals[i];

    /*  Free the variables you don't need anymore.                            */
    free(wrange_Index);
    free(wrange);
    free(Prange_Index);
    free(Prange);
    free(rho_legal);
    free(lambda_sky);

    if (dlp.start > dlp.arr_size){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rStarting index (start) is greater than the size of the array.\n"
        );
        return -1;
    }
    else if (dlp.start+dlp.n_used > dlp.arr_size){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rFinal index (start+n_used) is greater than size of array.\n"
        );
        return -1;
    }

    dlp.T_out = (complex double *)malloc((dlp.n_used+1)*sizeof(complex double));

    if (self->use_fft) DiffractionCorrectionSimpleFFT(&dlp);
    else {
        if (dlp.order == 0){
            if ((dlp.ecc == 0.0) && (dlp.peri == 0.0))
                if ((dlp.perturb[0] == 0) && (dlp.perturb[1] == 0) &&
                    (dlp.perturb[2] == 0) && (dlp.perturb[3] == 0) &&
                    (dlp.perturb[4] == 0)) {
                    DiffractionCorrectionNewton(&dlp);
                }
                else DiffractionCorrectionPerturbedNewton(&dlp);
            else DiffractionCorrectionEllipse(&dlp);
        }
        else if (dlp.order == 1) DiffractionCorrectionFresnel(&dlp);
        else DiffractionCorrectionLegendre(&dlp);
    }

    if (dlp.status == 0)
        dlp.n_used += 1;
    else if (dlp.status == 1){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rCould not extract data from inputs.\n"
        );
        return -1;
    }
    else if (dlp.status == 2){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\r\tRequired window width goes beyond the available data range.\n"
            "\r\t\tBad Point (Index):         \t%ld\n"
            "\r\t\tNumber of Points in Window:\t%ld\n"
            "\r\t\tDifference:                \t%ld\n"
            "\r\t\tSum:                       \t%ld\n"
            "\r\t\tArray Size:                \t%ld\n"
            "\r\tDifference must be positive and sum must\n"
            "\r\tbe less than array size.\n",
            dlp.start, dlp.n_used, dlp.start-dlp.n_used, dlp.start+dlp.n_used,
            dlp.arr_size
        );
        return -1;
    }
    else if (dlp.status == 3){
        PyErr_Format(
            PyExc_MemoryError,
            "\rError Encountered: rss_ringoccs"
            "\r\tspecial_functions.fresnel_transform\n\n"
            "\rMalloc failed to create new variables.\n"
            "\rYou are most likely out of memory.\n"
        );
        return -1;
    }
    else if (dlp.status == 4){
        PyErr_Format(
            PyExc_MemoryError,
            "\rError Encountered: rss_ringoccs"
            "\r\tspecial_functions.fresnel_transform\n\n"
            "\rInterp should be either 0, 2, 3, or 4."
        );
        return -1;
    }
    else if (dlp.status == 5){
        PyErr_Format(
            PyExc_MemoryError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\r\tIllegal string for wtype. Allowed string are:\n"
            "\r\t\trect\t\tRectangular Window\n"
            "\r\t\tcoss\t\tSquared Cosine Window\n"
            "\r\t\tkb20\t\tKaiser-Bessel with alpha=2.0 pi\n"
            "\r\t\tkb25\t\tKaiser-Bessel with alpha=2.5 pi\n"
            "\r\t\tkb35\t\tKaiser-Bessel with alpha=3.5 pi\n"
            "\r\t\tkbmd20\t\tModified Kaiser-Bessel with alpha=2.0 pi\n"
            "\r\t\tkbmd25\t\tModified Kaiser-Bessel with alpha=2.5 pi\n"
            "\r\t\tkbmd35\t\tModified Kaiser-Bessel with alpha=3.5 pi\n"
        );
        return -1;
    }
    else {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rInput should be one dimensional numpy arrays and\n"
            "\rreal valued numbers. The following are expected:\n"
            "\r\tT_in:        \t1-D complex numpy array\n"
            "\r\trho_km_vals: \t1-D real numpy array\n"
            "\r\tF_km_vals:   \t1-D real numpy array\n"
            "\r\tphi_rad_vals:\t1-D real numpy array\n"
            "\r\tkd_vals:     \t1-D real numpy array\n"
            "\r\tB_rad_vals:  \t1-D real numpy array\n"
            "\r\tD_km_vals:   \t1-D real numpy array\n"
            "\r\tw_km_vals:   \t1-D real numpy array\n"
            "\r\tstart:       \tPositive integer\n"
            "\r\tn_used:      \tPositive integer\n"
            "\r\twtype:       \tPositive integer\n"
            "\r\tuse_norm:    \tPositive integer\n"
            "\r\tuse_fwd:     \tPositive integer\n"
            "\r\torder:       \tPositive integer\n"
            "\r\tecc:         \tReal number\n"
            "\r\tperi:        \tReal number\n\n"
            "\rNOTE:\n"
            "\r\tOnly one dimensional numpy arrays are allowed. Only\n"
            "\r\tdouble types are supported. No current support for long\n"
            "\r\tdouble or float. Set this in Python with\n"
            "\r\tastype(numpy.float) or astype(numpy.float64).\n"
        );
        return -1;
    }
    RSS_RINGOCCSSetArray(self, dlp, rho_km_vals, tmp, NPY_DOUBLE, 
                         dlp.start, dlp.n_used);

    PyObject *T_out = PyArray_SimpleNewFromData(1, &dlp.n_used,
                                                NPY_CDOUBLE, (void *)dlp.T_out);

    tmp = self->T_vals;
    Py_INCREF(T_out);
    self->T_vals = T_out;
    Py_XDECREF(tmp);

    PyObject *T_in = PyArray_SimpleNewFromData(1, &dlp.n_used,
                                               NPY_CDOUBLE, (void *)dlp.T_in);

    tmp = self->T_hat_vals;
    Py_INCREF(T_in);
    self->T_hat_vals = T_in;
    Py_XDECREF(tmp);

    /*  Reset the phase to normal. This is very important! We modified the    *
     *  data in the dlp_inst directly, rather than steal a copy, so we need   *
     *  to put it back. Luckily, we only negated it. So just re-negate it.    */
    for (i=0; i<dlp.arr_size; ++i)
        dlp.phase_rad_vals[i] = -dlp.phase_rad_vals[i];

    /*  New free everything we've malloc'd, or we'll have memory leaks. Do    *
     *  NOT free the stuff malloc'd that we are storing in self. This will be *
     *  freed once the instance of the DiffractionCorrection class is         *
     *  destroyed. Only free the C-level stuff.                               */
    free(dlp.wtype);
    free(dlp.psitype);
    free(dlp.kd_vals);

    return 1;
}

RSS_RINGOCCSDiffModelingFourVars(gap_diffraction, Gap_Diffraction);
RSS_RINGOCCSDiffModelingFourVars(ringlet_diffraction, Ringlet_Diffraction);
RSS_RINGOCCSDiffModelingFourVarsREAL(ringlet_diffraction_phase,
                                     Ringlet_Diffraction_Phase);
RSS_RINGOCCSDiffModelingThreeVars(left_straightedge,
                                  Left_Straightedge_Diffraction);
RSS_RINGOCCSDiffModelingThreeVars(right_straightedge,
                                  Right_Straightedge_Diffraction);

static PyObject *square_wave_diffraction(PyObject *self, PyObject *args)
{
    /*  We'll need output and capsule for safely creating the output array    *
     *  and ensuring we don't have a memory leak. rho is the input variable.  */
    PyObject *output, *capsule, *rho;

    /*  Variables for the width of the waves and the Fresnel scale.           */
    double W, F;

    /*  The number of waves.                                                  */
    long N;

    /*  Variable for the size of the input array or list.                     */
    long dim;

    /*  Try to parse the user input, returning error if this fails.           */
    if (!PyArg_ParseTuple(args, "OddK", &rho, &W, &F, &N))
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.square_wave_diffraction\n\n"
            "\rCould not parse inputs. Legal inputs are:\n"
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"
            "\r\tW:     Positive constant (Float)\n"
            "\r\tF:     Positive constant (Float)\n"
            "\r\tN:     Positive Integer (Int)\n\n"
            "\rNotes:\n"
            "\r\trho must be a non-empty one dimensional numpy array."
        );
        return NULL;
    }

    /*  The width must be positive. Check this.                               */
    if (W <= 0.0)
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.square_wave_diffraction\n\n"
            "\rWidth is not positive. (i.e. W<=0)\n"
        );
        return NULL;
    }

    /*  As always, the Fresnel scale is a positive value. Check.              */
    else if (F <= 0.0)
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.square_wave_diffraction\n\n"
            "\rFresnel scale is not positive. (i.e. F<=0)\n"
        );
        return NULL;
    }

    /*  The number of waves to sum over must be positive as well.             */
    else if (N <= 0)
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.square_wave_diffraction\n\n"
            "\rThe number of waves is not positive (i.e. N<=0)\n"
        );
        return NULL;
    }

    /*  If the user supplied a float or int, simply pass the inputs to the    *
     *  C function and return the value back to Python.                       */
    if (PyLong_Check(rho) || PyFloat_Check(rho))
    {
        /*  Convert a python float or int into a C double.                    */
        double val = PyFloat_AsDouble(rho);

        /*  The output is a complex double, so declare one.                   */
        complex double out;

        /*  Now get the computed value and convert the C complex double to a  *
         *  Python complex value and return.                                  */
        out = Square_Wave_Diffraction_Double(val, W, F, N);

        /*  To build a Python complex value we need to pass the real and      *
         *  imaginary parts of out into PyComplex_FromDoubles. We can easily  *
         *  access these values using creal and cimag, which are functions    *
         *  found in complex.h.                                               */
        return PyComplex_FromDoubles(creal(out), cimag(out));                  
    }

    /*  As of v1.3 we allow lists to be passed to these routines.             */
    else if (PyList_Check(rho))
    {
        /*  If the user passed a list, we'll need to return one. Create a     *
         *  variable for indexing over the input list.                        */
        long i;

        /*  Create another variable for indexing over the list. This will be  *
         *  the object corresponding to the value of the index i.             */
        PyObject *ith_item;

        /*  And lastly, a double for storing the input value and a complex    *
         *  double for storing the output.                                    */
        double val;
        complex double out;

        /*  Get the number of elements in the list.                           */
        dim = PyList_Size(rho);

        /*  If the input list is empty, return with an error.                 */
        if (dim == 0)
        {
            PyErr_Format(
                PyExc_TypeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tdiffrec.square_wave_diffraction\n\n"
                "\rInput list is empty.\n"
            );
            return NULL;
        }

        /*  Create a new list to be returned to the user.                     */
        output = PyList_New(dim);

        /*  Loop over the entries of the list, examine values for error, and  *
         *  try to perform computations if legal inputs are given.            */
        for (i=0; i<dim; ++i)
        {
            /*  Extract the python object in the ith slot of the list.        */
            ith_item = PyList_GET_ITEM(rho, i);

            /*  The list should be homogeneous with real numbers only. Check  *
             *  each item and raise an error otherwise.                       */
            if (!PyFloat_Check(ith_item) && !PyLong_Check(ith_item))
            {
                PyErr_Format(PyExc_TypeError,
                             "\n\rError Encountered: rss_ringoccs\n"
                             "\r\tdiffrec.%s\n\n"
                             "\rInput list must contain real numbers only.\n");
                return NULL;
            }

            /*  Convert the python float to a C double.                       */
            val = PyFloat_AsDouble(ith_item);

            /*  Compute the value, set the item in the list, and move on to   *
             *  to the next one. We'll need to convert the complex double     *
             *  to a Python complex value as well. Again, we'll use creal and *
             *  cimag, found in complex.h, to perform this task.              */
            out = Square_Wave_Diffraction_Double(val, W, F, N);

            PyList_SET_ITEM(output, i,
                            PyComplex_FromDoubles(creal(out), cimag(out)));
        }

        /*  The list created without errors, so just return it.               */
        return output;
    }

    /*  If the input is not a float, int, or list, it must be a numpy array.  *
     *  Check this and return with error if not.                              */
    else if (!PyArray_Check(rho))
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.square_wave_diffraction\n\n"
            "\rCould not parse inputs. Legal inputs are:\n"
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"
            "\r\ta:     Positive constant (Float)\n"
            "\r\tb:     Positive constant (Float) greater than a\n"
            "\r\tF      Positive constant (Float)\n\n"
            "\rNotes:\n"
            "\r\trho must be a non-empty one dimensional numpy array,\n"
            "\r\ta float, or an int.\n"
        );
        return NULL;
    }

    /*  Variables needed for the array data.                                  */
    void *data;
    complex double *T_hat;

    /*  Try to convert rho to double. Raise an error otherwise.               */
    PyObject *newrho = PyArray_FromObject(rho, NPY_DOUBLE, 1, 1);

    /*  If PyArray_FromObject failed, newrho is NULL. Check this.             */
    if (!(newrho))
    {
        PyErr_Format(PyExc_TypeError,
                    "\n\rError Encountered: rss_ringoccs\n"
                    "\r\tdiffrec.square_wave_diffraction\n\n"
                    "\rInvalid data type for input array. Input should"
                    "\n\rbe a 1-dimensional array of real numbers.\n");
        return NULL;
    }

    /*  Check the inputs to make sure they're valid.                          */
    if (PyArray_NDIM((PyArrayObject *)newrho) != 1)
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.square_wave_diffraction\n\n"
            "\rInput numpy array is not one-dimensional.\n"
        );
        return NULL;
    }

    /*  Get the size of the input numpy array.                                */
    dim = PyArray_DIMS((PyArrayObject *)newrho)[0];

    /*  Check that the array isn't empty. Raise error otherwise.              */
    if (dim == 0)
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.square_wave_diffraction\n\n"
            "\rInput numpy array is empty.\n"
        );
        return NULL;
    }

    /*  Get a pointer to the actual data from the array. Allocate memory for  *
     *  the data of the output numpy array, which we'll call T_hat.           */
    data  = PyArray_DATA((PyArrayObject *)newrho);
    T_hat = (complex double *)malloc(dim*sizeof(complex double));

    /*  Variable for indexing over the array.                                 */
    long i;

    /*  Loop over the elements of the array and compute.                      */
    for (i=0; i<dim; ++i)
        T_hat[i] = Square_Wave_Diffraction_Double(((double *)data)[i], W, F, N);

    /*  Set the output and capsule, ensuring no memory leaks occur.           */
    output = PyArray_SimpleNewFromData(1, &dim, NPY_CDOUBLE, (void *)T_hat);
    capsule = PyCapsule_New((void *)T_hat, NULL, capsule_cleanup);

    /*  This frees the variable at the Python level once it's destroyed.      */
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);

    /*  Return the results to Python.                                         */
    return Py_BuildValue("N", output);
}

#define RSS_RINGOCCSCheckArrayAndSetPointer(arr, arr_value, arr_type, tmp, ptr,\
                                            dim, tempdim, incr)                \
/*  Check if the input is a python float or int. Store the data as a         */\
/*  double and set incr to 0 if it is.                                       */\
if (PyFloat_Check(arr) || PyLong_Check(arr))                                   \
{                                                                              \
    arr_value = PyFloat_AsDouble(arr);                                         \
    ptr = &arr_value;                                                          \
    incr = 0;                                                                  \
}                                                                              \
else                                                                           \
{                                                                              \
    tmp = PyArray_FromAny(arr, arr_type, 1, 1, NPY_ARRAY_BEHAVED, NULL);       \
                                                                               \
    if (!tmp)                                                                  \
    {                                                                          \
        PyErr_Format(PyExc_TypeError,                                          \
                    "\n\rError Encountered: rss_ringoccs\n"                    \
                    "\r\tdiffrec.fresnel_scale\n\n"                            \
                    "\rInvalid data type for one of the input arrays. Input\n" \
                    "\rshoule be a 1-dimensional array of real numbers.\n");   \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  Make sure the inputs are one dimensional objects.                    */\
    if (PyArray_NDIM((PyArrayObject *)tmp) != 1)                               \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.square_wave_diffraction\n\n"                          \
            "\rOne of the input numpy arrays is not one-dimensional.\n"        \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    tempdim = PyArray_DIMS((PyArrayObject *)tmp)[0];                           \
                                                                               \
    if (tempdim == 0)                                                          \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.square_wave_diffraction\n\n"                          \
            "\rOne of the input numpy arrays is empty.\n"                      \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
    else if (tempdim == 1)                                                     \
        incr = 0;                                                              \
    else                                                                       \
    {                                                                          \
        incr = 1;                                                              \
                                                                               \
        if (dim == 1)                                                          \
            dim = tempdim;                                                     \
        else if (dim != tempdim)                                               \
        {                                                                      \
            PyErr_Format(                                                      \
                PyExc_TypeError,                                               \
                "\n\rError Encountered: rss_ringoccs\n"                        \
                "\r\tdiffrec.square_wave_diffraction\n\n"                      \
                "\rTwo of the input arrays have different lengths.\n"          \
            );                                                                 \
            return NULL;                                                       \
        }                                                                      \
    }                                                                          \
                                                                               \
    ptr = (double *)PyArray_DATA((PyArrayObject *)tmp);                        \
}

static PyObject *fresnel_scale(PyObject *self, PyObject *args)
{
    /*  We'll need output and capsule for safely creating the output array    *
     *  and ensuring we don't have a memory leak. We'll also nneed PyObjects  *
     *  for each input variable.                                              */
    PyObject *output, *capsule, *lambda, *d_vals, *phi_vals, *b_vals;

    PyObject *tmp = Py_None;

    /*  Pointers for the data inside the arrays.                              */
    double *lambda_data, lambda_value;
    double *d_vals_data, d_value;
    double *phi_vals_data, phi_value;
    double *b_vals_data, b_value;

    /*  Description of the array type we need, which is numpy's double.       */
    PyArray_Descr *arr_type = PyArray_DescrFromType(NPY_DOUBLE);

    /*  Variables for the size of the input array.                            */
    long dim = 0;
    long tempdim = 0;

    /*  Variables for incrementing through the arrays.                        */
    unsigned char lambda_incr, d_incr, phi_incr, b_incr;

    /*  Try to parse the user input, returning error if this fails.           */
    if (!PyArg_ParseTuple(args, "OOOO", &lambda, &d_vals, &phi_vals, &b_vals))
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.fresnel_scale\n\n"
            "\rCould not parse inputs. Legal inputs are:\n"
            "\r\tlambda:   Numpy Array of positive real numbers (Floats)\n"
            "\r\tD:        Numpy Array of positive real numbers (Floats)\n"
            "\r\tphi:      Numpy Array of positive real numbers (Floats)\n"
            "\r\tB:        Numpy Array of positive real numbers (Floats)\n\n"
            "\rNotes:\n"
            "\r\tarrays must be a non-empty and one dimensional."
        );
        return NULL;
    }

    /*  Convert lambda to double, if possible.                                */
    if (PyFloat_Check(lambda) || PyLong_Check(lambda))
    {
        lambda_value = PyFloat_AsDouble(lambda);
        lambda_data  = &lambda_value;
        lambda_incr  = 0;
        dim = 1;
    }
    else
    {
        tmp = PyArray_FromAny(lambda, arr_type, 1, 1, NPY_ARRAY_BEHAVED, NULL);

        if (!tmp)
        {
            PyErr_Format(PyExc_TypeError,
                        "\n\rError Encountered: rss_ringoccs\n"
                        "\r\tdiffrec.fresnel_scale\n\n"
                        "\rInvalid data type for first input array. Input\n"
                        "\rshoule be a 1-dimensional array of real numbers.\n");
            return NULL;
        }

        /*  Make sure the inputs are one dimensional objects.                 */
        if (PyArray_NDIM((PyArrayObject *)tmp) != 1)
        {
            PyErr_Format(
                PyExc_TypeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tdiffrec.square_wave_diffraction\n\n"
                "\rFirst input numpy array is not one-dimensional.\n"
            );
            return NULL;
        }

        dim = PyArray_DIMS((PyArrayObject *)tmp)[0];

        if (dim == 0)
        {
            PyErr_Format(
                PyExc_TypeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tdiffrec.square_wave_diffraction\n\n"
                "\rFirst input numpy array is empty.\n"
            );
            return NULL;
        }

        lambda_data = (double *)PyArray_DATA((PyArrayObject *)tmp);
        lambda_incr = 1;
    }

    RSS_RINGOCCSCheckArrayAndSetPointer(d_vals, d_value, arr_type, tmp,
                                        d_vals_data, dim, tempdim, d_incr);

    RSS_RINGOCCSCheckArrayAndSetPointer(phi_vals, phi_value, arr_type, tmp,
                                        phi_vals_data, dim, tempdim, phi_incr);

    RSS_RINGOCCSCheckArrayAndSetPointer(b_vals, b_value, arr_type, tmp,
                                        b_vals_data, dim, tempdim, b_incr);

    /*  Allocate memory for the data of the output numpy array.               */
    double *F_vals_data = (double *)malloc(sizeof(double)*dim);

    /*  Variable for indexing over the array.                                 */
    long i;

    long j_lambda   = 0;
    long j_d_vals   = 0;
    long j_phi_vals = 0;
    long j_b_vals   = 0;

    /*  Loop over the elements of the array and compute.                      */
    for (i=0; i<dim; ++i)
    {
        F_vals_data[i] = Fresnel_Scale_Double(lambda_data[j_lambda],
                                              d_vals_data[j_d_vals],
                                              phi_vals_data[j_phi_vals],
                                              b_vals_data[j_b_vals]);

        j_lambda   += lambda_incr;
        j_d_vals   += d_incr;
        j_phi_vals += phi_incr;
        j_b_vals   += b_incr;
    }


    /*  Set the output and capsule, ensuring no memory leaks occur.           */
    output = PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE,
                                       (void *)F_vals_data);
    capsule = PyCapsule_New((void *)F_vals_data, NULL, capsule_cleanup);

    /*  This frees the variable at the Python level once it's destroyed.      */
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);

    /*  Return the results to Python.                                         */
    return Py_BuildValue("N", output);
}


static PyObject *
DiffractionCorrection(PyDiffrecObj *self, PyObject *Py_UNUSED(ignored))
{
    return PyUnicode_FromFormat("DiffractionCorrection");
}

static PyMethodDef diffrec_methods[] =
{
    {
        "square_wave_diffraction",
        square_wave_diffraction,
        METH_VARARGS,
        "Compute the normalized equivalent width of an array."
    },
    {
        "gap_diffraction",
        gap_diffraction,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "diffrec.gap_diffraction\n\r\t"
        "Purpose\n\r\t\t"
        "Compute the diffraction pattern of an annular gap in the plane.\n\r\t"
        "Arguments:\n\r\t\t"
        "rho (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of positive real numbers.\n\r\t\t"
        "a (float)\n\r\t\t\t"
        "A positive real number, the inner radius of the annulus.\n\r\t\t"
        "b (float)\n\r\t\t\t"
        "A positive real number, the outter radius of the annulus.\n\r\t\t"
        "F (float)\n\r\t\t\t"
        "A positive real number, the Fresnel scale (same units as rho).\n\r\t"
        "Outputs:\n\r\t\t"
        "T_hat (numpy.ndarray):\n\r\t\t\t"
        "Numpy array of complex numbers equal to the diffraction pattern.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import diffrec\n\r\t\t"
        ">>> x = numpy.arange(0,100,0.01)\n\r\t\t"
        ">>> a = 45\n\r\t\t"
        ">>> b = 55\n\r\t\t"
        ">>> F = 0.05\n\r\t\t"
        ">>> y = diffrec.gap_diffraction(x, 45, 55, 0.05)"
    },
    {
        "ringlet_diffraction",
        ringlet_diffraction,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "diffrec.gap_diffraction\n\r\t"
        "Purpose\n\r\t\t"
        "Compute the diffraction pattern of a ringlet in the plane.\n\r\t"
        "Arguments:\n\r\t\t"
        "rho (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of positive real numbers.\n\r\t\t"
        "a (float)\n\r\t\t\t"
        "A positive real number, the inner radius of the ring.\n\r\t\t"
        "b (float)\n\r\t\t\t"
        "A positive real number, the outter radius of the ring.\n\r\t\t"
        "F (float)\n\r\t\t\t"
        "A positive real number, the Fresnel scale (same units as rho).\n\r\t"
        "Outputs:\n\r\t\t"
        "T_hat (numpy.ndarray):\n\r\t\t\t"
        "Numpy array of complex numbers equal to the diffraction pattern.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import diffrec\n\r\t\t"
        ">>> x = numpy.arange(0,100,0.01)\n\r\t\t"
        ">>> a = 45\n\r\t\t"
        ">>> b = 55\n\r\t\t"
        ">>> F = 0.05\n\r\t\t"
        ">>> y = diffrec.ringlet_diffraction(x, a, b, F)"
    },
    {
        "ringlet_diffraction_phase",
        ringlet_diffraction_phase,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "diffrec.ringlet_diffraction_phase\n\r\t"
        "Purpose\n\r\t\t"
        "Compute the phase of a ringlet in the plane.\n\r\t"
        "Arguments:\n\r\t\t"
        "rho (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of positive real numbers.\n\r\t\t"
        "a (float)\n\r\t\t\t"
        "A positive real number, the inner radius of the ring.\n\r\t\t"
        "b (float)\n\r\t\t\t"
        "A positive real number, the outter radius of the ring.\n\r\t\t"
        "F (float)\n\r\t\t\t"
        "A positive real number, the Fresnel scale (same units as rho).\n\r\t"
        "Outputs:\n\r\t\t"
        "T_hat (numpy.ndarray):\n\r\t\t\t"
        "Numpy array of complex numbers equal to the diffraction pattern.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(0,100,0.01)\n\r\t\t"
        ">>> a = 45\n\r\t\t"
        ">>> b = 55\n\r\t\t"
        ">>> F = 0.05\n\r\t\t"
        ">>> y = diffrec.ringlet_diffraction_phase(x, a, b, F)"
    },
    {
        "right_straightedge",
        right_straightedge,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "diffrec.right_straightedge\n\r\t"
        "Purpose\n\r\t\t"
        "Compute the diffraction pattern of a straightedge.\n\r\t"
        "Arguments:\n\r\t\t"
        "rho (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of positive real numbers.\n\r\t\t"
        "a (float)\n\r\t\t\t"
        "A positive real number, starting point of the straightedge.\n\r\t\t"
        "F (float)\n\r\t\t\t"
        "A positive real number, the Fresnel scale (same units as rho).\n\r\t"
        "Outputs:\n\r\t\t"
        "T_hat (numpy.ndarray):\n\r\t\t\t"
        "Numpy array of complex numbers equal to the diffraction pattern.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import diffrec\n\r\t\t"
        ">>> x = numpy.arange(0,100,0.01)\n\r\t\t"
        ">>> a = 45\n\r\t\t"
        ">>> F = 0.05\n\r\t\t"
        ">>> y = diffrec.right_straightedge(x, a, F)"
    },
    {
        "left_straightedge",
        left_straightedge,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "diffrec.left_straightedge\n\r\t"
        "Purpose\n\r\t\t"
        "Compute the diffraction pattern of a straightedge.\n\r\t"
        "Arguments:\n\r\t\t"
        "rho (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of positive real numbers.\n\r\t\t"
        "a (float)\n\r\t\t\t"
        "A positive real number, starting point of the straightedge.\n\r\t\t"
        "F (float)\n\r\t\t\t"
        "A positive real number, the Fresnel scale (same units as rho).\n\r\t"
        "Outputs:\n\r\t\t"
        "T_hat (numpy.ndarray):\n\r\t\t\t"
        "Numpy array of complex numbers equal to the diffraction pattern.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import diffrec\n\r\t\t"
        ">>> x = numpy.arange(0,100,0.01)\n\r\t\t"
        ">>> a = 45\n\r\t\t"
        ">>> F = 0.05\n\r\t\t"
        ">>> y = diffrec.left_straightedge(x, a, F)"
    },
    {
        "fresnel_scale",
        fresnel_scale,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "diffrec.fresnel_scale\n\r\t"
        "Purpose\n\r\t\t"
        "Compute the Fresnel scale.\n\r\t"
        "Arguments:\n\r\t\t"
        "Lambda (numpy.ndarray):\n\r\t\t\t"
        "Numpy array of positive real numbers, the wavelength.\n\r\t\t"
        "D (numpy.ndarray)\n\r\t\t\t"
        "Numpy array of real numbers, distance to ring intercept point.\n\r\t\t"
        "phi (numpy.ndarray)\n\r\t\t\t"
        "Numpy array of real numbers, the ring azimuth angle (radians).\n\r\t"
        "B (numpy.ndarray)\n\r\t\t\t"
        "Numpy array of real numbers, the ring opening angle (radians).\n\r\t"
        "Outputs:\n\r\t\t"
        "F (numpy.ndarray):\n\r\t\t\t"
        "The Fresnel scale.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import diffrec\n\r\t\t"
        ">>> lambda = 3.6e-5"
        ">>> phi = numpy.arange(0,3.14,0.01)\n\r\t\t"
        ">>> B = 0.3\n\r\t\t"
        ">>> D = 2.0e6\n\r\t\t"
        ">>> y = diffrec.fresnel_scale(lambda, D, phi, B)"
    },
    {NULL, NULL, 0, NULL}
};

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

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "diffrec",
    NULL,
    -1,
    diffrec_methods,
    NULL,
    NULL,
    NULL,
    NULL
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