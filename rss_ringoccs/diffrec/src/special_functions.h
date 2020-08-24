#ifndef RSS_RINGOCCS_SPECIAL_FUNCTIONS_H
#define RSS_RINGOCCS_SPECIAL_FUNCTIONS_H

/*  complex data types, as well as _Complex_I, are defined here.              */
#include <complex.h>
#include <stdlib.h>

/*  Most functions return floats, not int or int-like data types. To allow    *
 *  most python/numpy data-types functions are created that intake the given  *
 *  data type, convert to double, and then return a double. This preprocessor *
 *  function eliminates excessive and repetitive code. This also ensures the  *
 *  naming convention of functions is consistent throughout the code. For     *
 *  example, consider the Sinc function defined by:                           *
 *                     --                                                     *
 *                    |     sin(x)/x ,      x =/= 0                           *
 *      sinc(x)  =    |                                                       *
 *                    |     1 ,             x = 0                             *
 *                     --                                                     *
 *  This routine is outlined in the sinc file where Sinc_Float, Sinc_Double,  *
 *  and Sinc_Long_Double are defined and use the C Standard Library functions *
 *  sinf, sin, and sinl which compute the sine function at float, double, and *
 *  long double precision accuracy, respectively. If a user at the Python     *
 *  level passes an array of integers it would be convenient to just have a   *
 *  function Sinc_Int that takes in a int, converts to double, and then       *
 *  returns as a double. Similarly if a user passes a char, long, short, etc. *
 *  This preprocessor routine generates all of these functions for a given    *
 *  function FuncName where FuncName_Double has ALREADY been defined in some  *
 *  external file, for example sinc.c.                                        */
#ifdef RSSRINGOCCSNonFloatInputForFloatOutput
#undef RSSRINGOCCSNonFloatInputForFloatOutput
#endif

#define RSSRINGOCCSNonFloatInputForFloatOutput(FuncName)\
double FuncName##_Char(char x)          {return FuncName##_Double((double)x);}\
double FuncName##_UChar(unsigned char x){return FuncName##_Double((double)x);}\
double FuncName##_Short(short x)        {return FuncName##_Double((double)x);}\
double FuncName##_UShort(unsigned short x){\
    return FuncName##_Double((double)x);\
}\
double FuncName##_Int(int x)            {return FuncName##_Double((double)x);}\
double FuncName##_UInt(unsigned int x)  {return FuncName##_Double((double)x);}\
double FuncName##_Long(long x)          {return FuncName##_Double((double)x);}\
double FuncName##_ULong(unsigned long x){return FuncName##_Double((double)x);}\
double FuncName##_Long_Long(long long x){return FuncName##_Double((double)x);}\
double FuncName##_ULong_Long(unsigned long long x){\
    return FuncName##_Double((double)x);\
}

/*  This code is for generating the code for the various window functions.    */
#ifdef RSSRINGOCCSNonFloatInputTwoVarForFloatOutput
#undef RSSRINGOCCSNonFloatInputTwoVarForFloatOutput
#endif

#define RSSRINGOCCSNonFloatInputTwoVarForFloatOutput(FuncName)\
double FuncName##_Char(char x, double y)                                    \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}                                                                           \
double FuncName##_UChar(unsigned char x, double y)                          \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}                                                                           \
double FuncName##_Short(short x, double y)                                  \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}                                                                           \
double FuncName##_UShort(unsigned short x, double y)                        \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}                                                                           \
double FuncName##_Int(int x, double y)                                      \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}                                                                           \
double FuncName##_UInt(unsigned int x, double y)                            \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}                                                                           \
double FuncName##_Long(long x, double y)                                    \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}                                                                           \
double FuncName##_ULong(unsigned long x, double y)                          \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}                                                                           \
double FuncName##_Long_Long(long long x, double y)                          \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}                                                                           \
double FuncName##_ULong_Long(unsigned long long x, double y)                \
{                                                                           \
    return FuncName##_Double((double)x, y);                                 \
}


/*  These lines are repeated over and over in this header file to define the  *
 *  various math functions. Use this preprocessor function to save lines of   *
 *  code, stay consistent with naming conventions, and add clarity.           */
#ifdef RSSRINGOCCSGenerateExternFunctions
#undef RSSRINGOCCSGenerateExternFunctions
#endif

#define RSSRINGOCCSGenerateExternFunctions(FuncName)                        \
extern float        FuncName##_Float(float x);                              \
extern double       FuncName##_Double(double x);                            \
extern long double  FuncName##_Long_Double(long double x);                  \
extern double       FuncName##_Char(char x);                                \
extern double       FuncName##_UChar(unsigned char x);                      \
extern double       FuncName##_Short(short x);                              \
extern double       FuncName##_UShort(unsigned short x);                    \
extern double       FuncName##_Int(int x);                                  \
extern double       FuncName##_UInt(unsigned int x);                        \
extern double       FuncName##_Long(long x);                                \
extern double       FuncName##_ULong(unsigned long x);                      \
extern double       FuncName##_Long_Long(long long x);                      \
extern double       FuncName##_ULong_Long(unsigned long long x);

/*  Macro to avoid repetition and ensure consistency in naming convention.    */
#ifdef RSSRINGOCCSTwoVarWindowFuncExtern
#undef RSSRINGOCCSTwoVarWindowFuncExtern
#endif

#define RSSRINGOCCSTwoVarWindowFuncExtern(FuncName)                         \
extern float        FuncName##_Float(float x, float W);                     \
extern double       FuncName##_Double(double x, double W);                  \
extern long double  FuncName##_Long_Double(long double x, long double W);   \
extern double       FuncName##_Char(char x, double W);                      \
extern double       FuncName##_UChar(unsigned char x, double W);            \
extern double       FuncName##_Short(short x, double W);                    \
extern double       FuncName##_UShort(unsigned short x, double W);          \
extern double       FuncName##_Int(int x, double W);                        \
extern double       FuncName##_UInt(unsigned int x, double W);              \
extern double       FuncName##_Long(long x, double W);                      \
extern double       FuncName##_ULong(unsigned long x, double W);            \
extern double       FuncName##_Long_Long(long long x, double W);            \
extern double       FuncName##_ULong_Long(unsigned long long x, double W);

#define RSSRINGOCCSThreeVarWindowFuncExtern(FuncName)                       \
extern float  FuncName##_Float(float x, float W, float alpha);              \
extern double FuncName##_Double(double x, double W, double alpha);          \
extern long double                                                          \
FuncName##_Long_Double(long double x, long double W, long double alpha);    \
extern double FuncName##_Char(char x, double W, double alpha);              \
extern double FuncName##_UChar(unsigned char x, double W, double alpha);    \
extern double FuncName##_Short(short x, double W, double alpha);            \
extern double FuncName##_UShort(unsigned short x, double W, double alpha);  \
extern double FuncName##_Int(int x, double W, double alpha);                \
extern double FuncName##_UInt(unsigned int x, double W, double alpha);      \
extern double FuncName##_Long(long x, double W, double alpha);              \
extern double FuncName##_ULong(unsigned long x, double W, double alpha);    \
extern double FuncName##_Long_Long(long long x, double W, double alpha);    \
extern double                                                               \
FuncName##_ULong_Long(unsigned long long x, double W, double alpha);

/*  Generate extern function names for all of the math functions.             */
RSSRINGOCCSGenerateExternFunctions(BesselJ0);
RSSRINGOCCSGenerateExternFunctions(BesselI0);
RSSRINGOCCSGenerateExternFunctions(LambertW);
RSSRINGOCCSGenerateExternFunctions(Sinc);
RSSRINGOCCSGenerateExternFunctions(Wavelength_To_Wavenumber);
RSSRINGOCCSGenerateExternFunctions(Frequency_To_Wavelength);
RSSRINGOCCSGenerateExternFunctions(Resolution_Inverse);
RSSRINGOCCSGenerateExternFunctions(Fresnel_Cos);
RSSRINGOCCSGenerateExternFunctions(Fresnel_Sin);

/******************************************************************************
 *------------------------Complex Fresnel Integral----------------------------*
 ******************************************************************************/
extern complex double Fresnel_Taylor_to_Asymptotic_Double(double x);
extern complex double Fresnel_While_to_Asymptotic_Func(double x);
extern complex double Fresnel_Heald_Rational_EPS_Minus_Three_Func(double x);
extern complex double Fresnel_Heald_Rational_EPS_Minus_Four_Func(double x);
extern complex double Fresnel_Heald_Rational_EPS_Minus_Six_Func(double x);
extern complex double Fresnel_Heald_Rational_EPS_Minus_Eight_Func(double x);

/*  Kaiser-Bessel function with alpha = 2pi, 2.5pi, and 3.5pi                 */
RSSRINGOCCSTwoVarWindowFuncExtern(Kaiser_Bessel_2_0);
RSSRINGOCCSTwoVarWindowFuncExtern(Kaiser_Bessel_2_5);
RSSRINGOCCSTwoVarWindowFuncExtern(Kaiser_Bessel_3_5);

/*  Modified Kaiser-Bessel with alpha = 2pi, 2.5pi, and 3.5pi                 */
RSSRINGOCCSTwoVarWindowFuncExtern(Modified_Kaiser_Bessel_2_0);
RSSRINGOCCSTwoVarWindowFuncExtern(Modified_Kaiser_Bessel_2_5);
RSSRINGOCCSTwoVarWindowFuncExtern(Modified_Kaiser_Bessel_3_5);

/*  Rectangular and Squared Cosine window functions.                          */
RSSRINGOCCSTwoVarWindowFuncExtern(Rect_Window);
RSSRINGOCCSTwoVarWindowFuncExtern(Coss_Window);

/*  Modified and unmodified Kaiser-Bessel function with arbitrary alpha.      */
RSSRINGOCCSThreeVarWindowFuncExtern(Kaiser_Bessel_Al);
RSSRINGOCCSThreeVarWindowFuncExtern(Modified_Kaiser_Bessel_Al);

#undef RSSRINGOCCSTwoVarWindowFuncExtern
#undef RSSRINGOCCSThreeVarWindowFuncExtern

extern void Legendre_Polynomials(double *legendre_p, double x, int order);

extern void
Alt_Legendre_Polynomials(double *poly, double *legendre_p, int order);

extern void
Fresnel_Kernel_Coefficients(double *fresnel_ker_coeffs, double *legendre_p,
                            double *alt_legendre_p, double Legendre_Coeff,
                            int order);

extern float        Max_Float(float *arr, long n_elements);
extern double       Max_Double(double *arr, long n_elements);
extern long         double Max_Long_Double(long double *arr, long n_elements);
extern short        Max_Short(short *arr, long n_elements);
extern int          Max_Int(int *arr, long n_elements);
extern long         Max_Long(long *arr, long n_elements);
extern long long    Max_Long_Long(long long *arr, long n_elements);

extern float        Min_Float(float *arr, long n_elements);
extern double       Min_Double(double *arr, long n_elements);
extern long         double Min_Long_Double(long double *arr, long n_elements);
extern short        Min_Short(short *arr, long n_elements);
extern int          Min_Int(int *arr, long n_elements);
extern long         Min_Long(long *arr, long n_elements);
extern long long    Min_Long_Long(long long *arr, long n_elements);

extern float        Normeq_Float(float *w_func, long n_elements);
extern double       Normeq_Double(double *w_func, long n_elements);
extern long double  Normeq_Long_Double(long double *w_func, long n_elements);
extern double       Normeq_Short(short *w_func, long n_elements);
extern double       Normeq_Int(int *w_func, long n_elements);
extern double       Normeq_Long(long *w_func, long n_elements);
extern double       Normeq_Long_Long(long long *w_func, long n_elements);

/*  Window Normalization Functions                                            */
extern float
Window_Normalization_Float(float *ker, long dim, float dx, float f_scale);

extern double
Window_Normalization_Double(double *ker, long dim, double dx, double f_scale);

extern long double
Window_Normalization_Long_Double(long double *ker, long dim,
                                 long double dx, long double f_scale);

extern double
Window_Normalization_Short(short *ker, long dim, double dx, double f_scale);

extern double
Window_Normalization_Int(int *ker, long dim, double dx, double f_scale);

extern double
Window_Normalization_Long(long *ker, long dim, double dx, double f_scale);

extern double
Window_Normalization_Long_Long(long long *ker, long dim,
                               double dx, double f_scale);

extern float
Window_Normalization_Complex_Float(complex float *ker, long dim,
                                   float dx, float f_scale);

extern double
Window_Normalization_Complex_Double(complex double *ker, long dim,
                                    double dx, double f_scale);

extern long double
Window_Normalization_Complex_Long_Double(complex long double *ker, long dim,
                                         long double dx, long double f_scale);

/*  Make sure these names are available.                                      */
#ifdef _get_one_real_from_three_real
#undef _get_one_real_from_three_real
#endif

#ifdef _get_complex_from_four_real
#undef _get_complex_from_four_real
#endif

#ifdef _get_complex_from_three_real
#undef _get_complex_from_three_real
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
#define _get_one_real_from_three_real(x1, x2, x3, y, dim, f) ({\
    /*  Declare necessary variables.                                         */\
    long i;\
    \
    for (i=0; i<dim; ++i){\
        y[i] = (*f)(x1[i], x2, x3);\
    }\
})

#define _get_complex_from_three_real(x, a, F, y, dim, f) ({\
    /*  Declare necessary variables.                                         */\
    long i;\
    \
    for (i=0; i<dim; ++i){\
        y[i] = (*f)(x[i], a, F);\
    }\
})

#define _get_complex_from_four_real(x, a, b, F, y, dim, f) ({\
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
        _get_one_real_from_one_real(((float *)data), y, dim, CName##_Float);\
        output  = PyArray_SimpleNewFromData(1, &dim, NPY_FLOAT, (void *)y);\
        capsule = PyCapsule_New(y, NULL, capsule_cleanup);\
    }\
    else if (typenum == NPY_DOUBLE){\
        double *y;\
        y = (double *)malloc(dim*sizeof(double));\
        _get_one_real_from_one_real(((double *)data), y, dim, CName##_Double);\
        output  = PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE, (void *)y);\
        capsule = PyCapsule_New(y, NULL, capsule_cleanup);\
    }\
    else if (typenum == NPY_LONGDOUBLE){\
        long double *y;\
        y = (long double *)malloc(dim*sizeof(long double));\
        _get_one_real_from_one_real(((long double *)data), y, dim,\
                                     CName##_Long_Double);\
        output = PyArray_SimpleNewFromData(1, &dim, NPY_LONGDOUBLE, (void *)y);\
        capsule = PyCapsule_New(y, NULL, capsule_cleanup);\
    }\
    else {\
        double *y;\
        y = (double *)malloc(dim*sizeof(double));\
        \
        if (typenum == NPY_BYTE)\
            _get_one_real_from_one_real(((char *)data), y, dim, CName##_Char);\
        else if (typenum == NPY_UBYTE)\
            _get_one_real_from_one_real(((unsigned char *)data),\
                                         y, dim, CName##_UChar);\
        else if (typenum == NPY_SHORT)\
            _get_one_real_from_one_real(((short *)data),\
                                         y, dim, CName##_Short);\
        else if (typenum == NPY_USHORT)\
            _get_one_real_from_one_real(((unsigned short *)data),\
                                         y, dim, CName##_UShort);\
        else if (typenum == NPY_INT)\
            _get_one_real_from_one_real(((int *)data), y, dim, CName##_Int);\
        else if (typenum == NPY_UINT)\
            _get_one_real_from_one_real(((unsigned int *)data),\
                                         y, dim, CName##_UInt);\
        else if (typenum == NPY_LONG)\
            _get_one_real_from_one_real(((long *)data), y, dim, CName##_Long);\
        else if (typenum == NPY_ULONG)\
            _get_one_real_from_one_real(((unsigned long *)data),\
                                         y, dim, CName##_ULong);\
        else if (typenum == NPY_LONGLONG)\
            _get_one_real_from_one_real(((long long *)data), y,\
                                         dim, CName##_Long_Long);\
        else if (typenum == NPY_ULONG)\
            _get_one_real_from_one_real(((unsigned long long *)data), y, dim,\
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
        _get_one_real_from_two_real(((float *)data), dx, y, dim,\
                                     CName##_Float);\
        output = PyArray_SimpleNewFromData(1, &dim, NPY_FLOAT, (void *)y);\
        capsule = PyCapsule_New(y, NULL, capsule_cleanup);\
    }\
    else if (typenum == NPY_LONGDOUBLE)\
    {\
        long double *y;\
        y = (long double *)malloc(dim*sizeof(long double));\
        _get_one_real_from_two_real(((long double *)data), dx, y, dim,\
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
            _get_one_real_from_two_real(((double *)data), dx, y, dim,\
                                         CName##_Double);\
        else if (typenum == NPY_BYTE)\
            _get_one_real_from_two_real(((char *)data), dx, y, dim,\
                                         CName##_Char);\
        else if (typenum == NPY_UBYTE)\
            _get_one_real_from_two_real(((unsigned char *)data), dx, y, dim,\
                                         CName##_UChar);\
        else if (typenum == NPY_SHORT)\
            _get_one_real_from_two_real(((short *)data), dx, y, dim,\
                                         CName##_Short);\
        else if (typenum == NPY_USHORT)\
            _get_one_real_from_two_real(((unsigned short *)data), dx, y, dim,\
                                         CName##_UShort);\
        else if (typenum == NPY_INT)\
            _get_one_real_from_two_real(((int *)data), dx, y, dim,\
                                         CName##_Int);\
        else if (typenum == NPY_UINT)\
            _get_one_real_from_two_real(((unsigned int *)data), dx, y, dim,\
                                         CName##_UInt);\
        else if (typenum == NPY_LONG)\
            _get_one_real_from_two_real(((long *)data), dx, y, dim,\
                                         CName##_Long);\
        else if (typenum == NPY_ULONG)\
            _get_one_real_from_two_real(((unsigned long *)data), dx, y, dim,\
                                         CName##_ULong);\
        else if (typenum == NPY_LONGLONG)\
            _get_one_real_from_two_real(((long long *)data), dx, y, dim,\
                                         CName##_Long_Long);\
        else if (typenum == NPY_ULONG)\
            _get_one_real_from_two_real(((unsigned long long *)data), dx, y,\
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

/*-------------Ringlet Diffraction Using Fresnel Approximation----------------*/

extern complex float
Ringlet_Diffraction_Float(float x, float a, float b, float F);

extern complex double
Ringlet_Diffraction_Double(double x, double a, double b, double F);

extern complex long double
Ringlet_Diffraction_Long_Double(long double x, long double a,
                                long double b, long double F);

extern double
Ringlet_Diffraction_Char(char x, double a, double b, double F);

extern double
Ringlet_Diffraction_UChar(unsigned char x, double a, double b, double F);

extern double
Ringlet_Diffraction_Short(short x, double a, double b, double F);

extern double
Ringlet_Diffraction_UShort(unsigned short x, double a, double b, double F);

extern double
Ringlet_Diffraction_Int(int x, double a, double b, double F);

extern double
Ringlet_Diffraction_UInt(unsigned int x, double a, double b, double F);

extern double
Ringlet_Diffraction_Long(long x, double a, double b, double F);

extern double
Ringlet_Diffraction_ULong(unsigned long x, double a, double b, double F);

extern double
Ringlet_Diffraction_Long_Long(long long x, double a, double b, double F);

extern double
Ringlet_Diffraction_ULong_Long(unsigned long long x, double a,
                               double b, double F);

/*----------------Gap Diffraction Using Fresnel Approximation-----------------*/

extern complex float
Gap_Diffraction_Float(float x, float a, float b, float F);

extern complex double
Gap_Diffraction_Double(double x, double a, double b, double F);

extern complex long double
Gap_Diffraction_Long_Double(long double x, long double a,
                            long double b, long double F);

extern double
Gap_Diffraction_Char(char x, double a, double b, double F);

extern double
Gap_Diffraction_UChar(unsigned char x, double a, double b, double F);

extern double
Gap_Diffraction_Short(short x, double a, double b, double F);

extern double
Gap_Diffraction_UShort(unsigned short x, double a, double b, double F);

extern double
Gap_Diffraction_Int(int x, double a, double b, double F);

extern double
Gap_Diffraction_UInt(unsigned int x, double a, double b, double F);

extern double
Gap_Diffraction_Long(long x, double a, double b, double F);

extern double
Gap_Diffraction_ULong(unsigned long x, double a, double b, double F);

extern double
Gap_Diffraction_Long_Long(long long x, double a, double b, double F);

extern double
Gap_Diffraction_ULong_Long(unsigned long long x, double a,
                           double b, double F);

/*-----------Ringlet Diffraction Phase Using Fresnel Approximation------------*/

extern float
Ringlet_Diffraction_Phase_Float(float x, float a, float b, float F);

extern double
Ringlet_Diffraction_Phase_Double(double x, double a, double b, double F);

extern long double
Ringlet_Diffraction_Phase_Long_Double(long double x, long double a,
                                      long double b, long double F);

/*--------Right Straight-Edge Diffraction Using Fresnel Approximation---------*/

extern complex float
Right_Straightedge_Diffraction_Float(float x, float edge, float F);

extern complex double
Right_Straightedge_Diffraction_Double(double x, double edge, double F);

extern complex long double
Right_Straightedge_Diffraction_Long_Double(long double x, long double edge,
                                           long double F);

extern double
Right_Straightedge_Diffraction_Char(char x, double a, double F);

extern double
Right_Straightedge_Diffraction_UChar(unsigned char x, double a, double F);

extern double
Right_Straightedge_Diffraction_Short(short x, double a, double F);

extern double
Right_Straightedge_Diffraction_UShort(unsigned short x, double a, double F);

extern double
Right_Straightedge_Diffraction_Int(int x, double a, double F);

extern double
Right_Straightedge_Diffraction_UInt(unsigned int x, double a, double F);

extern double
Right_Straightedge_Diffraction_Long(long x, double a, double F);

extern double
Right_Straightedge_Diffraction_ULong(unsigned long x, double a, double F);

extern double
Right_Straightedge_Diffraction_Long_Long(long long x, double a, double F);

extern double
Right_Straightedge_Diffraction_ULong_Long(unsigned long long x,
                                          double a, double F);

/*---------Left Straight-Edge Diffraction Using Fresnel Approximation---------*/

extern complex float
Left_Straightedge_Diffraction_Float(float x, float edge, float F);

extern complex double
Left_Straightedge_Diffraction_Double(double x, double edge, double F);

extern complex long double
Left_Straightedge_Diffraction_Long_Double(long double x, long double edge,
                                          long double F);

extern complex float
Left_Straightedge_Diffraction_Float(float x, float edge, float F);

extern complex double
Left_Straightedge_Diffraction_Double(double x, double edge, double F);

extern complex long double
Left_Straightedge_Diffraction_Long_Double(long double x, long double edge,
                                          long double F);

extern double
Left_Straightedge_Diffraction_Char(char x, double a, double F);

extern double
Left_Straightedge_Diffraction_UChar(unsigned char x, double a, double F);

extern double
Left_Straightedge_Diffraction_Short(short x, double a, double F);

extern double
Left_Straightedge_Diffraction_UShort(unsigned short x, double a, double F);

extern double
Left_Straightedge_Diffraction_Int(int x, double a, double F);

extern double
Left_Straightedge_Diffraction_UInt(unsigned int x, double a, double F);

extern double
Left_Straightedge_Diffraction_Long(long x, double a, double F);

extern double
Left_Straightedge_Diffraction_ULong(unsigned long x, double a, double F);

extern double
Left_Straightedge_Diffraction_Long_Long(long long x, double a, double F);

extern double
Left_Straightedge_Diffraction_ULong_Long(unsigned long long x,
                                         double a, double F);

/*------------Square Wave Diffraction Using Fresnel Approximation-------------*/
extern complex double
Square_Wave_Diffraction_Double(double x, double W, double F, long N);

#endif
