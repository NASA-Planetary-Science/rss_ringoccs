#ifndef RSS_RINGOCCS_SPECIAL_FUNCTIONS_H
#define RSS_RINGOCCS_SPECIAL_FUNCTIONS_H

/*  complex data types, as well as _Complex_I, are defined here.              */
#include <complex.h>
#include <stdlib.h>

/*  compute_norm_eq, max and min found here. math.h included here as well.    */
#include "__fresnel_diffraction.h"
#include "__math_functions.h"
#include "__get_array.h"
#include "__window_functions.h"
#include "__diffraction_functions.h"

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

#ifdef OneVarFunctionForNumpy
#undef OneVarFunctionForNumpy
#endif

#ifdef WindowFunctionForNumpy
#undef WindowFunctionForNumpy
#endif

#ifdef VarToString
#undef VarToString
#endif

#define VarToString(Var) (#Var)

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

/*  Again, to avoid repetition, the code used in the API for numpy is arrays  *
 *  is identical for all functions with the exception of the function name.   *
 *  This preprocessor function saves thousands of lines of code, all of which *
 *  would be just copy/paste otherwise.                                       */
#define OneVarFunctionForNumpy(FuncName, CName)\
static PyObject * FuncName(PyObject *self, PyObject *args)\
{\
    /*  Declare necessary variables.                                         */\
    PyObject *output, *capsule;\
    PyArrayObject *x;\
    char typenum;\
    long dim;\
    void *data;\
    \
    /*  Parse the data from Python and try to convert it to a usable format. */\
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &x)){\
        PyErr_Format(PyExc_TypeError,\
                     "\n\rError Encountered: rss_ringoccs\n"\
                     "\r\tdiffrec.special_functions.%s\n\n"\
                     "\rCould not parse inputs. Legal inputs are:\n"\
                     "\r\tx: Numpy Array of real numbers (Floats)\n\rNotes:\n"\
                     "\r\tx must be a non-empty one dimensional numpy array.",\
                     VarToString(FuncName));\
        return NULL;\
    }\
    \
    /*  Grab useful information about the data.                              */\
    typenum = (char)PyArray_TYPE(x);\
    dim     = PyArray_DIMS(x)[0];\
    data    = PyArray_DATA(x);\
    \
    /*  Check the inputs to make sure they're valid.                         */\
    if (PyArray_NDIM(x) != 1){\
        PyErr_Format(PyExc_TypeError, "\n\rError Encountered: rss_ringoccs\n"\
                                      "\r\tdiffrec.special_functions.%s\n"\
                                      "\n\rInput is not 1-dimensional.\n",\
                                      VarToString(FuncName));\
        return NULL;\
    }\
    else if (dim == 0){\
        PyErr_Format(PyExc_TypeError, "\n\rError Encountered: rss_ringoccs\n"\
                                      "\r\tdiffrec.special_functions.%s"\
                                      "\n\n\rInput numpy array is empty.\n",\
                                      VarToString(FuncName));\
    }\
    \
    if (typenum == NPY_FLOAT){\
        float *y;\
        y = (float *)malloc(dim*sizeof(float));\
        __get_one_real_from_one_real(((float *)data), y, dim, CName##_Float);\
        output  = PyArray_SimpleNewFromData(1, &dim, NPY_FLOAT, (void *)y);\
        capsule = PyCapsule_New(y, NULL, capsule_cleanup);\
    }\
    else if (typenum == NPY_DOUBLE){\
        double *y;\
        y = (double *)malloc(dim*sizeof(double));\
        __get_one_real_from_one_real(((double *)data), y, dim, CName##_Double);\
        output  = PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE, (void *)y);\
        capsule = PyCapsule_New(y, NULL, capsule_cleanup);\
    }\
    else if (typenum == NPY_LONGDOUBLE){\
        long double *y;\
        y = (long double *)malloc(dim*sizeof(long double));\
        __get_one_real_from_one_real(((long double *)data), y, dim,\
                                     CName##_Long_Double);\
        output = PyArray_SimpleNewFromData(1, &dim, NPY_LONGDOUBLE, (void *)y);\
        capsule = PyCapsule_New(y, NULL, capsule_cleanup);\
    }\
    else {\
        double *y;\
        y = (double *)malloc(dim*sizeof(double));\
        \
        if (typenum == NPY_BYTE)\
            __get_one_real_from_one_real(((char *)data), y, dim, CName##_Char);\
        else if (typenum == NPY_UBYTE)\
            __get_one_real_from_one_real(((unsigned char *)data),\
                                         y, dim, CName##_UChar);\
        else if (typenum == NPY_SHORT)\
            __get_one_real_from_one_real(((short *)data),\
                                         y, dim, CName##_Short);\
        else if (typenum == NPY_USHORT)\
            __get_one_real_from_one_real(((unsigned short *)data),\
                                         y, dim, CName##_UShort);\
        else if (typenum == NPY_INT)\
            __get_one_real_from_one_real(((int *)data), y, dim, CName##_Int);\
        else if (typenum == NPY_UINT)\
            __get_one_real_from_one_real(((unsigned int *)data),\
                                         y, dim, CName##_UInt);\
        else if (typenum == NPY_LONG)\
            __get_one_real_from_one_real(((long *)data), y, dim, CName##_Long);\
        else if (typenum == NPY_ULONG)\
            __get_one_real_from_one_real(((unsigned long *)data),\
                                         y, dim, CName##_ULong);\
        else if (typenum == NPY_LONGLONG)\
            __get_one_real_from_one_real(((long long *)data), y,\
                                         dim, CName##_Long_Long);\
        else if (typenum == NPY_ULONG)\
            __get_one_real_from_one_real(((unsigned long long *)data), y, dim,\
                                         CName##_Long_Long);\
        else {\
            PyErr_Format(PyExc_TypeError,\
                         "\n\rError Encountered: rss_ringoccs\n"\
                         "\r\tdiffrec.special_functions.%s\n\n"\
                         "\rInvalid data type for input array. Input should be"\
                         "\n\ra 1-dimensional array of real numbers.\n",\
                         VarToString(FuncName));\
            return NULL;\
        }\
        \
        output  = PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE, (void *)y);\
        capsule = PyCapsule_New(y, NULL, capsule_cleanup);\
    }\
    \
    /*  This frees the variable at the Python level once it's destroyed.     */\
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);\
    \
    /*  Return the results to Python.                                        */\
    return Py_BuildValue("N", output);\
}

#define WindowFunctionForNumpy(FuncName, CName)\
static PyObject * FuncName(PyObject *self, PyObject *args)\
{\
    PyObject *output, *capsule;\
    PyArrayObject *x;\
    double dx;\
    char typenum;\
    long dim;\
    void *data;\
    \
    if (!PyArg_ParseTuple(args, "O!d", &PyArray_Type, &x, &dx))\
    {\
        PyErr_Format(PyExc_TypeError,\
                     "\n\rError Encountered: rss_ringoccs\n"\
                     "\r\tdiffrec.special_functions.%s\n\n"\
                     "\rCould not parse inputs. Legal inputs are:\n"\
                     "\r\tx:     Numpy Array of real numbers (Floats)\n"\
                     "\r\tdx:    Positive real number (Float)\n\rNotes:\n"\
                     "\r\tx must be a non-empty one dimensional numpy array.",\
                     VarToString(FuncName));\
        return NULL;\
    }\
    \
    /*  Useful information about the data.                                   */\
    typenum = (char)PyArray_TYPE(x);\
    dim     = PyArray_DIMS(x)[0];\
    data    = PyArray_DATA(x);\
    \
    /*  Check the inputs to make sure they're valid.                         */\
    if (PyArray_NDIM(x) != 1)\
    {\
        PyErr_Format(PyExc_TypeError, "\n\rError Encountered: rss_ringoccs\n"\
                                      "\r\tdiffrec.special_functions.%s\n\n"\
                                      "\rInput array is not 1-dimensional.\n",\
                                      VarToString(FuncName));\
        return NULL;\
    }\
    else if (dim == 0)\
    {\
        PyErr_Format(PyExc_TypeError, "\n\rError Encountered: rss_ringoccs\n"\
                                      "\r\tdiffrec.special_functions.%s\n\n"\
                                      "\rInput numpy array is empty.\n",\
                                      VarToString(FuncName));\
        return NULL;\
    }\
    \
    /*  Check that dx is positive.                                           */\
    if (dx <= 0)\
    {\
        PyErr_Format(PyExc_ValueError, "\n\rError Encountered: rss_ringoccs\n"\
                                       "\r\tdiffrec.special_functions.%s\n\n"\
                                       "\rdx must be a positive number.\n",\
                                       VarToString(FuncName));\
        return NULL;\
    }\
    \
    if (typenum == NPY_FLOAT)\
    {\
        float *y;\
        y = (float *)malloc(dim*sizeof(float));\
        __get_one_real_from_two_real(((float *)data), dx, y, dim,\
                                     CName##_Float);\
        output = PyArray_SimpleNewFromData(1, &dim, NPY_FLOAT, (void *)y);\
        capsule = PyCapsule_New(y, NULL, capsule_cleanup);\
    }\
    else if (typenum == NPY_LONGDOUBLE)\
    {\
        long double *y;\
        y = (long double *)malloc(dim*sizeof(long double));\
        __get_one_real_from_two_real(((long double *)data), dx, y, dim,\
                                     CName##_Long_Double);\
        output = PyArray_SimpleNewFromData(1, &dim, NPY_LONGDOUBLE, (void *)y);\
        capsule = PyCapsule_New(y, NULL, capsule_cleanup);\
    }\
    else\
    {\
        double *y;\
        y = (double *)malloc(dim*sizeof(double));\
        \
        if (typenum == NPY_DOUBLE)\
            __get_one_real_from_two_real(((double *)data), dx, y, dim,\
                                         CName##_Double);\
        else if (typenum == NPY_BYTE)\
            __get_one_real_from_two_real(((char *)data), dx, y, dim,\
                                         CName##_Char);\
        else if (typenum == NPY_UBYTE)\
            __get_one_real_from_two_real(((unsigned char *)data), dx, y, dim,\
                                         CName##_UChar);\
        else if (typenum == NPY_SHORT)\
            __get_one_real_from_two_real(((short *)data), dx, y, dim,\
                                         CName##_Short);\
        else if (typenum == NPY_USHORT)\
            __get_one_real_from_two_real(((unsigned short *)data), dx, y, dim,\
                                         CName##_UShort);\
        else if (typenum == NPY_INT)\
            __get_one_real_from_two_real(((int *)data), dx, y, dim,\
                                         CName##_Int);\
        else if (typenum == NPY_UINT)\
            __get_one_real_from_two_real(((unsigned int *)data), dx, y, dim,\
                                         CName##_UInt);\
        else if (typenum == NPY_LONG)\
            __get_one_real_from_two_real(((long *)data), dx, y, dim,\
                                         CName##_Long);\
        else if (typenum == NPY_ULONG)\
            __get_one_real_from_two_real(((unsigned long *)data), dx, y, dim,\
                                         CName##_ULong);\
        else if (typenum == NPY_LONGLONG)\
            __get_one_real_from_two_real(((long long *)data), dx, y, dim,\
                                         CName##_Long_Long);\
        else if (typenum == NPY_ULONG)\
            __get_one_real_from_two_real(((unsigned long long *)data), dx, y,\
                                         dim, CName##_ULong_Long);\
        else\
        {\
            PyErr_Format(PyExc_TypeError,\
                         "\n\rError Encountered: rss_ringoccs\n"\
                         "\r\tdiffrec.special_functions.%s\n\n"\
                         "\rInvalid data type for input array. Input should\n"\
                         "\rbe a 1-dimensional numpy array of real numbers.\n",\
                         VarToString(FuncName));\
            return NULL;\
        }\
        \
        output = PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE, (void *)y);\
        capsule = PyCapsule_New(y, NULL, capsule_cleanup);\
    }\
    \
    /*  This frees the variable at the Python level once it's destroyed.     */\
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);\
    \
    /*  Return the results to Python.                                        */\
    return Py_BuildValue("N", output);\
}

/******************************************************************************
 *--------------------Single Slit Fraunhofer Diffraction----------------------*
 ******************************************************************************/

extern float
Single_Slit_Fraunhofer_Diffraction_Float(float x, float z, float a);

extern double
Single_Slit_Fraunhofer_Diffraction_Double(double x, double z, double a);

extern long double
Single_Slit_Fraunhofer_Diffraction_Long_Double(long double x, long double z,
                                               long double a);

/******************************************************************************
 *--------------------Double Slit Fraunhofer Diffraction----------------------*
 ******************************************************************************/

extern float
Double_Slit_Fraunhofer_Diffraction_Float(float x, float z, float a,
                                         float d, float lambda);

extern double
Double_Slit_Fraunhofer_Diffraction_Double(double x, double z, double a,
                                          double d, double lambda);

extern long double
Double_Slit_Fraunhofer_Diffraction_Long_Double(long double x, long double z,
                                               long double a, long double d,
                                               long double lambda);

/******************************************************************************
 *------------------------------Fresnel Scale---------------------------------*
 ******************************************************************************/

extern float
Fresnel_Scale_Float(float lambda, float d, float phi, float b);

extern double
Fresnel_Scale_Double(double lambda, double d, double phi, double b);

extern long double
Fresnel_Scale_Long_Double(long double lambda, long double d,
                          long double phi, long double b);

/******************************************************************************
 *--------------------------Frequency to Wavelength---------------------------*
 ******************************************************************************/

extern float        Frequency_To_Wavelength_Float(float frequency);
extern double       Frequency_To_Wavelength_Double(double frequency);
extern long double  Frequency_To_Wavelength_Long_Double(long double frequency);

/******************************************************************************
 *-------------------------Wavelength to Wavenumber---------------------------*
 ******************************************************************************/

extern float        Wavelength_To_Wavenumber_Float(float frequency);
extern double       Wavelength_To_Wavenumber_Double(double frequency);
extern long double  Wavelength_To_Wavenumber_Long_Double(long double frequency);

/******************************************************************************
 *-----------------------------------Where------------------------------------*
 ******************************************************************************/

extern long **
Where_Greater_Char(char *data, long dim, double threshold);

extern long **
Where_Greater_UChar(unsigned char *data, long dim, double threshold);

extern long **
Where_Greater_Short(short *data, long dim, double threshold);

extern long **
Where_Greater_UShort(unsigned short *data, long dim, double threshold);

extern long **
Where_Greater_Int(int *data, long dim, double threshold);

extern long **
Where_Greater_UInt(unsigned int *data, long dim, double threshold);

extern long **
Where_Greater_Long(long *data, long dim, double threshold);

extern long **
Where_Greater_ULong(unsigned long *data, long dim, double threshold);

extern long **
Where_Greater_Long_Long(long long *data, long dim, double threshold);

extern long **
Where_Greater_ULong_Long(unsigned long long *data, long dim, double threshold);

extern long **
Where_Greater_Float(float *data, long dim, double threshold);

extern long **
Where_Greater_Double(double *data, long dim, double threshold);

extern long **
Where_Greater_Long_Double(long double *data, long dim, long double threshold);

extern long **
Where_Lesser_Char(char *data, long dim, double threshold);

extern long **
Where_Lesser_UChar(unsigned char *data, long dim, double threshold);

extern long **
Where_Lesser_Short(short *data, long dim, double threshold);

extern long **
Where_Lesser_UShort(unsigned short *data, long dim, double threshold);

extern long **
Where_Lesser_Int(int *data, long dim, double threshold);

extern long **
Where_Lesser_UInt(unsigned int *data, long dim, double threshold);

extern long **
Where_Lesser_Long(long *data, long dim, double threshold);

extern long **
Where_Lesser_ULong(unsigned long *data, long dim, double threshold);

extern long **
Where_Lesser_Long_Long(long long *data, long dim, double threshold);

extern long **
Where_Lesser_ULong_Long(unsigned long long *data, long dim, double threshold);

extern long **
Where_Lesser_Float(float *data, long dim, double threshold);

extern long **
Where_Lesser_Double(double *data, long dim, double threshold);

extern long **
Where_Lesser_Long_Double(long double *data, long dim, long double threshold);

#endif