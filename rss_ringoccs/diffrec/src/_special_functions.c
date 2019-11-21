/*  To avoid compiler warnings about deprecated numpy stuff.                  */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/*  complex data types, as well as _Complex_I, are defined here.              */
#include <complex.h>

/*  compute_norm_eq, max and min found here. math.h included here as well.    */
#include "__math_functions.h"
#include "__where.h"
#include "__get_array.h"

/*  Window functions and Fresnel transforms defined here.                     */
#include "__diffraction_functions.h"

/* All wrapper functions defined within these files.                          */
#include "_fraunhofer_diffraction_wrappers.h"
#include "__fresnel_diffraction.h"
#include "_fresnel_integrals_wrappers.h"
#include "_fresnel_kernel_wrappers.h"
#include "_lambertw_wrappers.h"
#include "_physics_functions_wrappers.h"
#include "_resolution_inverse_function_wrappers.h"
#include "_sinc_wrappers.h"

/*  Various header files required for the C-Python API to work.               */
#include <Python.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

/*  Make sure the name __get_where_pointer is available.                      */
#ifdef __get_one_real_from_one_real
#undef __get_one_real_from_one_real
#endif

#ifdef __get_complex_from_four_real
#undef __get_complex_from_four_real
#endif

#ifdef __get_complex_from_three_real
#undef __get_complex_from_three_real
#endif

/*  To avoid repeating the same code over and over again, define this macro   *
 *  to be used for all of the where_lesser functions. Since the only thing    *
 *  that changes between the various functions is the type of the input       *
 *  pointer, the code is exactly the same.                                    */

#define __get_one_real_from_one_real(x, y, dim, f) ({\
    /*  Declare necessary variables.                                         */\
    long i;\
    \
    for (i=0; i<dim; ++i){\
        y[i] = (*f)(x[i]);\
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

/*---------------------------DEFINE PYTHON FUNCTIONS--------------------------*
 *  This contains the Numpy-C and Python-C API parts that allow for the above *
 *  functions to be called in Python. Numpy arrays, as well as floating point *
 *  and integer valued arguments may then be passed into these functions for  *
 *  improvement in performance, as opposed to the routines written purely in  *
 *  Python. Successful compiling requires the Numpy and Python header files.  *
 *----------------------------------------------------------------------------*/

/*  This function frees the memory allocated to a pointer by malloc when the  *
 *  corresponding variable is destroyed at the Python level.                  */
static void capsule_cleanup(PyObject *capsule)
{
    void *memory = PyCapsule_GetPointer(capsule, NULL);
    free(memory);
}

/******************************************************************************
 *  Function:                                                                 *
 *      where_greater                                                         *
 *  Purpose:                                                                  *
 *      Given a numpy array arr, and a real number threshold, returns an      *
 *      array of alll indices n such that arr[n] > threshold.                 *
 *  Arguments:                                                                *
 *      arr (Numpy Array):                                                    *
 *          An arbitrary numpy array of real numbers.                         *
 *      threshold (double):                                                   *
 *          A real valued number. This is the value such that, for all n such *
 *          that arr[n] > threshold, the index n will be appended to the      *
 *          returning array.                                                  *
 *  Output:                                                                   *
 *      where_arr:                                                            *
 *          A numpy array of integer corresponding to the indices n such that *
 *          arr[n] > threshold.                                               *
 *  NOTES:                                                                    *
 *      1.) Values that are equal to threshold will not be included.          *
 *      2.) The input array MUST be one-dimensional.                          *
 *      3.) integer and floating point types are allowed, but complex are not.*
 ******************************************************************************/
static PyObject *where_greater(PyObject *self, PyObject *args)
{
    /*  Declare necessary variables.                                          */
    PyObject *output, *capsule;
    PyArrayObject *arr;
    double threshold;
    long **where;

    /*  The input is a numpy array, proceed. Otherwise spit out an error.     */
    if (PyArg_ParseTuple(args, "O!d", &PyArray_Type, &arr, &threshold)){
        
        /*  Declare some more necessary variables.                            */
        npy_int typenum, dim;
        void *data;

        /*  Check to make sure input is one dimensional.                      */
        if (PyArray_NDIM(arr) != 1){
            PyErr_Format(
                PyExc_TypeError,
                "rss_ringoccs.diffrec.special_functions.where_greater\n"
                "\r\tInput must be a one-dimensional array and a real number."
            );
            return NULL;
        }

        /*   Useful information about the data.                               */
        typenum = PyArray_TYPE(arr);
        dim     = PyArray_DIMS(arr)[0];
        data    = PyArray_DATA(arr);

        /*  Compute the index array corresponding to the indices n such that  *
         *  arr[n] > threshold. The returned pointer is created using malloc  *
         *  so we must remember to free it later to avoid memory leaks.       */
        switch(typenum)
        {
            case NPY_BYTE:
                where = Where_Greater_Char((char *)data, dim, threshold);
                break;
            case NPY_UBYTE:
                where =
                Where_Greater_UChar((unsigned char *)data, dim, threshold);
                break;
            case NPY_SHORT:
                where = Where_Greater_Short((short *)data, dim, threshold);
                break;
            case NPY_USHORT:
                where =
                Where_Greater_UShort((unsigned short *)data, dim, threshold);
                break;
            case NPY_INT:
                where = Where_Greater_Int((int *)data, dim, threshold);
                break;
            case NPY_UINT:
                where =
                Where_Greater_UInt((unsigned int *)data, dim, threshold);
                break;
            case NPY_LONG:
                where = Where_Greater_Long((long *)data, dim, threshold);
                break;
            case NPY_ULONG:
                where =
                Where_Greater_ULong((unsigned long *)data, dim, threshold);
                break;
            case NPY_LONGLONG:
                where =
                Where_Greater_Long_Long((long long *)data, dim, threshold);
                break;
            case NPY_ULONGLONG:
                where = Where_Greater_ULong_Long((unsigned long long *)data,
                                                 dim, threshold);
                break;
            case NPY_FLOAT:
                where = Where_Greater_Float((float *)data, dim, threshold);
                break;
            case NPY_DOUBLE:
                where = Where_Greater_Double((double *)data, dim, threshold);
                break;
            case NPY_LONGDOUBLE:
                where =
                Where_Greater_Long_Double((long double *)data, dim, threshold);
                break;
            default:
                PyErr_Format(
                    PyExc_TypeError,
                    "\n\r\trss_ringoccs.diffrec.math_functions.where_greater\n"
                    "\r\t\tInput numpy array should be real valued."
                );
                return NULL;
        }

        /*  The first element of where is the array of indices.               */
        long *where_arr = where[0];

        /*  The second element is the number of points in the index array.    */
        long where_dim = *where[1];

        /*  Create a Numpy array object to be passed back to Python.          */
        output =
        PyArray_SimpleNewFromData(1, &where_dim, NPY_LONG, (void *)where_arr);

        /*  This frees the variable at the Python level once it's destroyed.  */
        capsule = PyCapsule_New(where_arr, NULL, capsule_cleanup);
        PyArray_SetBaseObject((PyArrayObject *)output, capsule);

        /*  Where is created using malloc, so make sure to free it.           */
        free(where);

        /*  Return the results to Python.                                     */
        return Py_BuildValue("N", output);
    }
    else {
        /*  If he input is not a numpy array, return to Python with an error. */
        PyErr_Format(
            PyExc_TypeError,
            "\n\r\trss_ringoccs.diffrec.special_functions.where_greater\n"
            "\r\t\tInput should be a real numpy array and a real number."
        );
        return NULL;
    }
}

/******************************************************************************
 *  Function:                                                                 *
 *      where_lesser                                                          *
 *  Purpose:                                                                  *
 *      Given a numpy array arr, and a real number threshold, returns an      *
 *      array of all indices n such that arr[n] < threshold.                  *
 *  Arguments:                                                                *
 *      arr (Numpy Array):                                                    *
 *          An arbitrary numpy array of real numbers.                         *
 *      threshold (double):                                                   *
 *          A real valued number. This is the value such that, for all n such *
 *          that arr[n] < threshold, the index n will be appended to the      *
 *          returning array.                                                  *
 *  Output:                                                                   *
 *      where_arr:                                                            *
 *          A numpy array of integer corresponding to the indices n such that *
 *          arr[n] < threshold.                                               *
 *  NOTES:                                                                    *
 *      1.) Values that are equal to threshold will not be included.          *
 *      2.) The input array MUST be one-dimensional.                          *
 *      3.) integer and floating point types are allowed, but complex are not.*
 ******************************************************************************/
static PyObject *where_lesser(PyObject *self, PyObject *args)
{
    /*  Declare necessary variables.                                          */
    PyObject *output, *capsule;
    PyArrayObject *arr;
    double threshold;
    long **where;

    /*  The input is a numpy array, proceed. Otherwise spit out an error.     */
    if (PyArg_ParseTuple(args, "O!d", &PyArray_Type, &arr, &threshold)){

        /*  Check to make sure input is one dimensional.                      */
        if (PyArray_NDIM(arr) != 1){
            PyErr_Format(
                PyExc_TypeError,
                "rss_ringoccs.diffrec.special_functions.where_greater\n"
                "\r\tInput must be a one-dimensional array and a real number."
            );
            return NULL;
        }

        /*   Useful information about the data.                               */
        npy_int typenum = PyArray_TYPE(arr);
        npy_int dim     = PyArray_DIMS(arr)[0];

        /*  Compute the index array corresponding to the indices n such that  *
         *  arr[n] > threshold. The returned pointer is created using malloc  *
         *  so we must remember to free it later to avoid memory leaks.       */
        if (typenum == NPY_BYTE){
            char *data = (char *)PyArray_DATA(arr);
            where = Where_Lesser_Char(data, dim, threshold);
        }
        else if (typenum == NPY_UBYTE){
            unsigned char *data = (unsigned char *)PyArray_DATA(arr);
            where = Where_Lesser_UChar(data, dim, threshold);
        }
        else if (typenum == NPY_SHORT){
            short *data = (short *)PyArray_DATA(arr);
            where = Where_Lesser_Short(data, dim, threshold);
        }
        else if (typenum == NPY_USHORT){
            unsigned short *data = (unsigned short *)PyArray_DATA(arr);
            where = Where_Lesser_UShort(data, dim, threshold);
        }
        else if (typenum == NPY_INT){
            int *data = (int *)PyArray_DATA(arr);
            where = Where_Lesser_Int(data, dim, threshold);
        }
        else if (typenum == NPY_UINT){
            unsigned int *data = (unsigned int *)PyArray_DATA(arr);
            where = Where_Lesser_UInt(data, dim, threshold);
        }
        else if (typenum == NPY_LONG){
            long *data = (long *)PyArray_DATA(arr);
            where = Where_Lesser_Long(data, dim, threshold);
        }
        else if (typenum == NPY_ULONG){
            unsigned long *data = (unsigned long *)PyArray_DATA(arr);
            where = Where_Lesser_ULong(data, dim, threshold);
        }
        else if (typenum == NPY_LONGLONG){
            long long *data = (long long *)PyArray_DATA(arr);
            where = Where_Lesser_Long_Long(data, dim, threshold);
        }
        else if (typenum == NPY_ULONG){
            unsigned long long *data = (unsigned long long *)PyArray_DATA(arr);
            where = Where_Lesser_ULong_Long(data, dim, threshold);
        }
        else if (typenum == NPY_FLOAT){
            float *data = (float *)PyArray_DATA(arr);
            where = Where_Lesser_Float(data, dim, threshold);
        }
        else if (typenum == NPY_DOUBLE){
            double *data = (double *)PyArray_DATA(arr);
            where = Where_Lesser_Double(data, dim, threshold);
        }
        else if (typenum == NPY_LONGDOUBLE){
            long double *data = (long double *)PyArray_DATA(arr);
            where = Where_Lesser_Long_Double(data, dim, threshold);
        }
        else {
            /*  If he input is not a numpy array, return to Python with an error. */
            PyErr_Format(
                PyExc_TypeError,
                "\n\r\trss_ringoccs.diffrec.math_functions.where_greater\n"
                "\r\t\tInput numpy array should be real valued."
            );
            return NULL;
        }

        /*  The first element of where is the array of indices.               */
        long *where_arr = where[0];

        /*  The second element is the number of points in the index array.    */
        long where_dim = *where[1];

        /*  Create a Numpy array object to be passed back to Python.          */
        output  = PyArray_SimpleNewFromData(1, &where_dim, NPY_LONG,
                                            (void *)where_arr);

        /*  This frees the variable at the Python level once it's destroyed.  */
        capsule = PyCapsule_New(where_arr, NULL, capsule_cleanup);
        PyArray_SetBaseObject((PyArrayObject *)output, capsule);

        /*  Where is created using malloc, so make sure to free it.           */
        free(where);

        /*  Return the results to Python.                                     */
        return Py_BuildValue("N", output);
    }
    else {
        /*  If he input is not a numpy array, return to Python with an error. */
        PyErr_Format(
            PyExc_TypeError,
            "\n\r\trss_ringoccs.diffrec.math_functions.where_greater\n"
            "\r\t\tInput should be a numpy array of numbers and a real number."
        );
        return NULL;
    }
}

static PyObject *compute_norm_eq(PyObject *self, PyObject *args)
{
    PyArrayObject *arr;
    PyObject *tuple = PyTuple_GetItem(args, 0);

    if (PyLong_Check(tuple)){
        long normeq;
        if (PyArg_ParseTuple(args, "l", &normeq)){
            return PyLong_FromLong(normeq);
        }
        else {
            PyErr_Format(
                PyExc_TypeError,
                "\n\r\trss_ringoccs.diffrec.math_functions.compute_norm_eq\n"
                "\r\t\tCould not parse int type input.");
            return NULL;
        }
    }
    else if (PyFloat_Check(tuple)){
        double normeq;
        if (PyArg_ParseTuple(args, "d", &normeq)){
            return PyFloat_FromDouble(normeq);
        }
        else {
            PyErr_Format(
                PyExc_TypeError,
                "\n\r\trss_ringoccs.diffrec.math_functions.compute_norm_eq\n"
                "\r\t\tCould not parse float type input.");
            return NULL;
        }
    }
    else if (PyArg_ParseTuple(args, "O!", &PyArray_Type, &arr)){
        npy_int typenum, dim;
        void *data;

        // Check to make sure input isn't zero dimensional!
        if (PyArray_NDIM(arr) != 1){
            PyErr_Format(
                PyExc_TypeError,
                "\n\trss_ringoccs.diffrec.special_functions.compute_norm_eq\n"
                "\r\t\tInput must be a one-dimensional array."
            );
            return NULL;
        }

        // Useful information about the data.
        typenum = PyArray_TYPE(arr);
        dim     = PyArray_DIMS(arr)[0];
        data    = PyArray_DATA(arr);
        if (dim == 0){
            PyErr_Format(
                PyExc_TypeError,
                "\n\r\trss_ringoccs.diffrec.math_functions.compute_norm_eq\n"
                "\r\t\tInput is zero dimensional."
            );
            return NULL;
        }
        if (typenum == NPY_FLOAT){
            return PyFloat_FromDouble(
                Normeq_Float((float *)data, dim)
            );
        }
        else if (typenum == NPY_DOUBLE){
            return PyFloat_FromDouble(
                Normeq_Double((double *)data, dim)
            );
        }
        else if (typenum == NPY_LONGDOUBLE){
            return PyFloat_FromDouble(
                Normeq_Long_Double((long double *)data, dim)
            );
        }
        else if (typenum == NPY_SHORT){
            return PyFloat_FromDouble(
                Normeq_Short((short *)data, dim)
            );
        }
        else if (typenum == NPY_INT){
            return PyFloat_FromDouble(
                Normeq_Int((int *)data, dim)
            );
        }
        else if (typenum == NPY_LONG){
            return PyFloat_FromDouble(
                Normeq_Long((long *)data, dim)
            );
        }
        else if (typenum == NPY_LONGLONG){
            return PyFloat_FromDouble(
                Normeq_Long_Long((long long *)data, dim)
            );
        }
        else {
            PyErr_Format(
                PyExc_TypeError,
                "\n\r\trss_ringoccs.diffrec.math_functions.compute_norm_eq\n"
                "\r\t\tInput should be a numpy array of numbers."
            );
            return NULL;
        }
    }
    else {
        PyErr_Format(
            PyExc_TypeError,
            "\n\r\trss_ringoccs.diffrec.math_functions.compute_norm_eq\n"
            "\r\t\tInput should be a numpy array of numbers."
        );
        return NULL;
    }
}

static PyObject *max(PyObject *self, PyObject *args)
{
    PyArrayObject *arr;
    PyObject *tuple = PyTuple_GetItem(args, 0);

    if (PyLong_Check(tuple)){
        long max;
        if (PyArg_ParseTuple(args, "l", &max)){
            return PyLong_FromLong(max);
        }
        else {
            PyErr_Format(PyExc_TypeError,
                         "\n\r\trss_ringoccs.diffrec.special_functions.max\n"
                         "\r\t\tCould not parse int type input.");
            return NULL;
        }
    }
    else if (PyFloat_Check(tuple)){
        double max;
        if(PyArg_ParseTuple(args, "d", &max)){
            return PyFloat_FromDouble(max);
        }
        else {
            PyErr_Format(PyExc_TypeError,
                         "\n\r\trss_ringoccs.diffrec.special_functions.max\n"
                         "\r\t\tCould not parse float type input.");
            return NULL;
        }
    }
    else if (PyArg_ParseTuple(args, "O!", &PyArray_Type, &arr)){
        npy_int typenum, dim;
        void *data;

        // Check to make sure input isn't zero dimensional!
        if (PyArray_NDIM(arr) != 1) goto FAIL;

        // Useful information about the data.
        typenum = PyArray_TYPE(arr);
        dim     = PyArray_DIMS(arr)[0];
        data    = PyArray_DATA(arr);

        if (dim == 0) goto FAIL;

        switch(typenum)
        {
            case NPY_FLOAT:
                return PyFloat_FromDouble(Max_Float((float *)data, dim));
            case NPY_DOUBLE:
                return PyFloat_FromDouble(Max_Double((double *)data, dim));
            case NPY_LONGDOUBLE:
                return
                PyFloat_FromDouble(Max_Long_Double((long double *)data, dim));
            case NPY_SHORT:
                return PyLong_FromLong(Max_Short((short *)data, dim));
            case NPY_INT:
                return PyLong_FromLong(Max_Int((int *)data, dim));
            case NPY_LONG:
                return PyLong_FromLong(Max_Long((long *)data, dim));
            case NPY_LONGLONG:
                return PyLong_FromLong(Max_Long_Long((long long *)data, dim));
            default:
                goto FAIL;
        }
    }
    else goto FAIL;
    FAIL: {
        PyErr_Format(PyExc_TypeError,
                     "\n\r\trss_ringoccs.diffrec.special_functions.max\n"
                     "\r\t\tInput should be a one dimensional numpy array of\n"
                     "\r\t\treal numbers, or a float/int number.\n"
                     "\r\t\tExample:\n"
                     "\r\t\t\t>>> import numpy\n"
                     "\r\t\t\t>>> import _special_functions as sf\n"
                     "\r\t\t\t>>> x = numpy.random.rand(100)\n"
                     "\r\t\t\t>>> y = sf.max(x)\n\n"
                     "\r\t\tNOTE:\n"
                     "\r\t\t\tOnly one dimensional numpy arrays are allowed.\n"
                     "\r\t\t\tComplex numbers are not allowed. If the input\n"
                     "\r\t\t\tis a single floating point or integer number,\n"
                     "\r\t\t\tthe output will simply be that number.");
        return NULL;
    };
}

static PyObject *min(PyObject *self, PyObject *args)
{
    PyArrayObject *arr;
    PyObject *tuple = PyTuple_GetItem(args, 0);

    if (PyLong_Check(tuple)){
        long min;
        if (PyArg_ParseTuple(args, "l", &min)){
            return PyLong_FromLong(min);
        }
        else {
            PyErr_Format(PyExc_TypeError,
                         "\n\r\trss_ringoccs.diffrec.math_functions.min\n"
                         "\r\t\tCould not parse int type input.");
            return NULL;
        }
    }
    else if (PyFloat_Check(tuple)){
        double min;
        if(PyArg_ParseTuple(args, "d", &min)){
            return PyFloat_FromDouble(min);
        }
        else {
            PyErr_Format(PyExc_TypeError,
                         "\n\r\trss_ringoccs.diffrec.math_functions.min\n"
                         "\r\t\tCould not parse float type input.");
            return NULL;
        }
    }
    else if (PyArg_ParseTuple(args, "O!", &PyArray_Type, &arr)){
        npy_int typenum, dim;
        void *data;

        // Check to make sure input isn't zero dimensional!
        if (PyArray_NDIM(arr) != 1){
            PyErr_Format(PyExc_TypeError,
                         "\n\r\trss_ringoccs.diffrec.math_functions.min\n"
                         "\r\t\tInput must be one dimensional.");
            return NULL;
        }

        // Useful information about the data.
        typenum = PyArray_TYPE(arr);
        dim     = PyArray_DIMS(arr)[0];
        data    = PyArray_DATA(arr);

        if (dim == 0){
            PyErr_Format(PyExc_TypeError,
                         "\n\r\trss_ringoccs.diffrec.math_functions.min\n"
                         "\r\t\tInput is zero dimensional.");
            return NULL;
        }

        if (typenum == NPY_FLOAT){
            return PyFloat_FromDouble(
                Min_Float((float *)data, dim)
            );
        }
        else if (typenum == NPY_DOUBLE){
            return PyFloat_FromDouble(
                Min_Double((double *)data, dim)
            );
        }
        else if (typenum == NPY_LONGDOUBLE){
            return PyFloat_FromDouble(
                Min_Long_Double((long double *)data, dim)
            );
        }
        else if (typenum == NPY_SHORT){
            return PyLong_FromLong(
                Min_Short((short *)data, dim)
            );
        }
        else if (typenum == NPY_INT){
            return PyLong_FromLong(
                Min_Int((int *)data, dim)
            );
        }
        else if (typenum == NPY_LONG){
            return PyLong_FromLong(
                Min_Long((long *)data, dim)
            );
        }
        else if (typenum == NPY_LONGLONG){
            return PyLong_FromLong(
                Min_Long_Long((long long *)data, dim)
            );
        }
        else {
            PyErr_Format(PyExc_TypeError,
                         "\n\r\trss_ringoccs.diffrec.math_functions.min\n"
                         "\r\t\tInput should be a numpy array of numbers.");
            return NULL;
        }
    }
    else {
        PyErr_Format(PyExc_TypeError,
                     "rss_ringoccs.diffrec.math_functions.min\n"
                     "\n\r\trss_ringoccs.diffrec.math_functions.min\n"
                     "\r\t\tInput should be a numpy array of numbers.");
        return NULL;
    }
}

static PyObject *fresnel_transform(PyObject *self, PyObject *args)
{
    PyObject *output, *capsule;
    npy_int dim;
    PyArrayObject *T_in,    *rho_km_vals, *F_km_vals, *phi_rad_vals;
    PyArrayObject *kd_vals, *B_rad_vals,  *D_km_vals, *w_km_vals;
    long start, n_used;
    unsigned char wtype, use_norm, use_fwd, order;
    double ecc, peri;


    if (!(PyArg_ParseTuple(args, "O!O!O!O!O!O!O!O!llbbbbdd",
                           &PyArray_Type, &T_in,
                           &PyArray_Type, &rho_km_vals,
                           &PyArray_Type, &F_km_vals,
                           &PyArray_Type, &phi_rad_vals,
                           &PyArray_Type, &kd_vals,
                           &PyArray_Type, &B_rad_vals,
                           &PyArray_Type, &D_km_vals,
                           &PyArray_Type, &w_km_vals,
                           &start, &n_used, &wtype, &use_norm,
                           &use_fwd, &order, &ecc, &peri))){ 
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
        return NULL;
    }

    DLPObj dlp;

    dlp.start    = start;
    dlp.n_used   = n_used;
    dlp.wtype    = wtype;
    dlp.use_norm = use_norm;
    dlp.use_fwd  = use_fwd;
    dlp.order    = order;
    dlp.ecc      = ecc;
    dlp.peri     = peri;


    dim = PyArray_DIMS(T_in)[0];

    /*  Error checks on T_in.                                                 */
    if (PyArray_TYPE(T_in) != NPY_CDOUBLE){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rFirst argument (T_in) must be a complex numpy array.\n"
        );
        return NULL;
    }
    else if (PyArray_NDIM(T_in) != 1){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rFirst argument (T_in) must be one dimensional.\n"
            "\r\tNumber of dimensions: %ld", PyArray_NDIM(T_in)
        );
        return NULL;
    }
    else dlp.T_in = (complex double *)PyArray_DATA(T_in);

    /*  Error checks on rho_km_vals                                           */
    if (PyArray_TYPE(rho_km_vals) != NPY_DOUBLE){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rSecond argument (rho_km_vals) must be a real numpy array.\n"
        );
        return NULL;
    }
    else if (PyArray_NDIM(rho_km_vals) != 1){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rSecond argument (rho_km_vals) must be one dimensional.\n"
        );
        return NULL;
    }
    else if (PyArray_DIMS(rho_km_vals)[0] != dim){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rrho_km_vals and T_in have a different number of elements.\n"
        );
        return NULL;
    }
    else dlp.rho_km_vals = (double *)PyArray_DATA(rho_km_vals);

    /*  Error checks on F_km_vals.                                            */
    if (PyArray_TYPE(F_km_vals) != NPY_DOUBLE){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rThird argument (F_km_vals) must be a real numpy array.\n"
        );
        return NULL;
    }
    else if (PyArray_NDIM(F_km_vals) != 1){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rThird argument (F_km_vals) must be one dimensional.\n"
        );
        return NULL;
    }
    else if (PyArray_DIMS(F_km_vals)[0] != dim){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rF_km_vals and T_in have a different number of elements.\n"
        );
        return NULL;
    }
    else dlp.F_km_vals = (double *)PyArray_DATA(F_km_vals);

    /*  Error checks on phi_rad_vals.                                         */
    if (PyArray_TYPE(phi_rad_vals) != NPY_DOUBLE){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rFourth argument (phi_rad_vals) must be a real numpy array.\n"
        );
        return NULL;
    }
    else if (PyArray_NDIM(phi_rad_vals) != 1){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rFourth argument (phi_rad_vals) must be one dimensional.\n"
        );
        return NULL;
    }
    else if (PyArray_DIMS(phi_rad_vals)[0] != dim){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rphi_rad_vals and T_in have a different number of elements.\n"
        );
        return NULL;
    }
    else dlp.phi_rad_vals = (double *)PyArray_DATA(phi_rad_vals);

    /*  Error checks on kd_vals.                                              */
    if (PyArray_TYPE(kd_vals) != NPY_DOUBLE){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rFifth argument (kd_vals) must be a real numpy array.\n"
        );
        return NULL;
    }
    else if (PyArray_NDIM(kd_vals) != 1){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rFifth argument (kd_vals) must be one dimensional.\n"
        );
        return NULL;
    }
    else if (PyArray_DIMS(kd_vals)[0] != dim){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rkd_vals and T_in have a different number of elements.\n"
        );
        return NULL;
    }
    else dlp.kd_vals = (double *)PyArray_DATA(kd_vals);

    /*  Error checks on B_rad_vals.                                           */
    if (PyArray_TYPE(B_rad_vals) != NPY_DOUBLE){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rSixth argument (B_rad_vals) must be a real numpy array.\n"
        );
        return NULL;
    }
    else if (PyArray_NDIM(B_rad_vals) != 1){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rSixth argument (B_rad_vals) must be one dimensional.\n"
        );
        return NULL;
    }
    else if (PyArray_DIMS(B_rad_vals)[0] != dim){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rB_rad_vals and T_in have a different number of elements.\n"
        );
        return NULL;
    }
    else dlp.B_rad_vals = (double *)PyArray_DATA(B_rad_vals);

    /*  Error checks on D_km_vals.                                            */
    if (PyArray_TYPE(D_km_vals) != NPY_DOUBLE){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rSeventh argument (D_km_vals) must be a real numpy array.\n"
        );
        return NULL;
    }
    else if (PyArray_NDIM(D_km_vals) != 1){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rSeventh argument (D_km_vals) must be one dimensional.\n"
        );
        return NULL;
    }
    else if (PyArray_DIMS(B_rad_vals)[0] != dim){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rD_km_vals and T_in have a different number of elements.\n"
        );
        return NULL;
    }
    else dlp.D_km_vals = (double *)PyArray_DATA(D_km_vals);

    /*  Error checks on w_km_vals.                                            */
    if (PyArray_TYPE(w_km_vals) != NPY_DOUBLE){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rEigth argument (w_km_vals) must be a real numpy array.\n"
        );
        return NULL;
    }
    else if (PyArray_NDIM(w_km_vals) != 1){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rEigth argument (w_km_vals) must be one dimensional.\n"
        );
        return NULL;
    }
    else if (PyArray_DIMS(w_km_vals)[0] != dim){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rw_km_vals and T_in have a different number of elements.\n"
        );
        return NULL;
    }
    else dlp.w_km_vals = (double *)PyArray_DATA(w_km_vals);

    if (start > dim){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rStarting index (start) is greater than the size of the array.\n"
        );
        return NULL;
    }
    else if (start+n_used > dim){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rFinal index (start+n_used) is greater than size of array.\n"
        );
        return NULL;
    }

    dlp.T_out = (complex double *)malloc((dlp.n_used+1)*sizeof(complex double));

    if (dlp.order == 0){
        if ((dlp.ecc == 0.0) && (dlp.peri == 0.0))
            DiffractionCorrectionNewton(&dlp);
        else DiffractionCorrectionEllipse(&dlp);
    }
    else if (dlp.order == 1) DiffractionCorrectionFresnel(&dlp);
    else DiffractionCorrectionLegendre(&dlp);

    if (dlp.status == 0){

        /*  Create a Numpy array object to be passed back to Python.          */
        dlp.n_used += 1;
        output = PyArray_SimpleNewFromData(1, &dlp.n_used,
                                           NPY_CDOUBLE, (void *)dlp.T_out);

        /*  This frees the variable at the Python level once it's destroyed.  */
        capsule = PyCapsule_New(dlp.T_out, NULL, capsule_cleanup);
        PyArray_SetBaseObject((PyArrayObject *)output, capsule);

        /*  Return the results to Python.                                     */
        return Py_BuildValue("N", output);
    }
    else if (dlp.status == 1){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\rCould not extract data from inputs.\n"
        );
        return NULL;
    }
    else if (dlp.status == 2){
        PyErr_Format(
            PyExc_IndexError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.fresnel_transform\n\n"
            "\r\tRequired window width goes beyond the available data range.\n"
            "\r\t\tStarting Point:            \t%ld\n"
            "\r\t\tNumber of Points in Window:\t%ld\n"
            "\r\t\tDifference:                \t%ld\n\n"
            "\r\tDifference must be positive.\n",
            dlp.start, dlp.n_used, dlp.start-dlp.n_used
        );
        return NULL;
    }
    else if (dlp.status == 3){
        PyErr_Format(
            PyExc_MemoryError,
            "\rError Encountered: rss_ringoccs"
            "\r\tspecial_functions.fresnel_transform\n\n"
            "\rMalloc failed to create new variables.\n"
            "\rYou are most likely out of memory.\n"
        );
        return NULL;
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
        return NULL;
    }
}

static PyObject *gap_diffraction(PyObject *self, PyObject *args)
{
    PyObject *output, *capsule;
    PyArrayObject *rho;
    double a, b, F;
    char typenum;
    long dim;
    void *data;

    if (!PyArg_ParseTuple(args, "O!ddd", &PyArray_Type, &rho, &a, &b, &F)){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.gap_diffraction\n\n"
            "\rCould not parse inputs. Legal inputs are:\n"
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"
            "\r\ta:     Positive constant (Float)\n"
            "\r\tb:     Positive constant (Float) greater than a\n"
            "\r\tF      Positive constant (Float)\n\n"
            "\rNotes:\n"
            "\r\trho must be a non-empty one dimensional numpy array."
        );
        return NULL;
    }

    /*  Useful information about the data.                                */
    typenum = (char)PyArray_TYPE(rho);
    dim     = PyArray_DIMS(rho)[0];
    data    = PyArray_DATA(rho);

    /*  Check the inputs to make sure they're valid.                          */
    if (PyArray_NDIM(rho) != 1){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.gap_diffraction\n\n"
            "\rInput numpy array is not one-dimensional.\n"
        );
        return NULL;
    }
    else if (a >= b){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.gap_diffraction\n\n"
            "\rInner radius is not less than outer radius (i.e. a >= b).\n"
        );
        return NULL;
    }
    else if (a <= 0.0){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.gap_diffraction\n\n"
            "\rInner radius is negative. (i.e. a<0)\n"
        );
        return NULL;
    }
    else if (F < 0.0){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.gap_diffraction\n\n"
            "\rFresnel scale is negative (i.e. F<0).\n"
        );
        return NULL;
    }
    else if (dim == 0){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.gap_diffraction\n\n"
            "\rInput numpy array is empty.\n"
        );
    }

    if (typenum == NPY_FLOAT){
        complex float *T_hat;
        T_hat = (complex float *)malloc(dim*sizeof(complex float));
        __get_complex_from_four_real(((float *)data), a, b, F, T_hat,
                                     dim, Gap_Diffraction_Float);
        output = PyArray_SimpleNewFromData(1, &dim, NPY_CFLOAT, (void *)T_hat);
        capsule = PyCapsule_New(T_hat, NULL, capsule_cleanup);
    }
    else if (typenum == NPY_DOUBLE){
        complex double *T_hat;
        T_hat = (complex double *)malloc(dim*sizeof(complex double));
        __get_complex_from_four_real(((double *)data), a, b, F, T_hat, dim,
                                     Gap_Diffraction_Double);
        output = PyArray_SimpleNewFromData(1, &dim, NPY_CDOUBLE, (void *)T_hat);
        capsule = PyCapsule_New(T_hat, NULL, capsule_cleanup);
    }
    else if (typenum == NPY_LONGDOUBLE){
        complex long double *T_hat;
        T_hat = (complex long double *)malloc(dim*sizeof(complex long double));
        __get_complex_from_four_real(((long double *)data), a, b, F, T_hat,
                                     dim, Gap_Diffraction_Long_Double);
        output = PyArray_SimpleNewFromData(1, &dim, NPY_CLONGDOUBLE,
                                           (void *)T_hat);
        capsule = PyCapsule_New(T_hat, NULL, capsule_cleanup);
    }
    else {
        complex double *T_hat;
        T_hat = (complex double *)malloc(dim*sizeof(complex double));

        if (typenum == NPY_BYTE){
            __get_complex_from_four_real(((char *)data), a, b, F, T_hat, dim,
                                         Gap_Diffraction_Char);
        }
        else if (typenum == NPY_UBYTE){
            __get_complex_from_four_real(((unsigned char *)data), a, b, F,
                                         T_hat, dim, Gap_Diffraction_UChar);
        }
        else if (typenum == NPY_SHORT){
            __get_complex_from_four_real(((short *)data), a, b, F, T_hat, dim,
                                         Gap_Diffraction_Short);
        }
        else if (typenum == NPY_USHORT){
            __get_complex_from_four_real(((unsigned short *)data), a, b, F,
                                         T_hat, dim, Gap_Diffraction_UShort);
        }
        else if (typenum == NPY_INT){
            __get_complex_from_four_real(((int *)data), a, b, F, T_hat, dim,
                                         Gap_Diffraction_Int);
        }
        else if (typenum == NPY_UINT){
            __get_complex_from_four_real(((unsigned int *)data), a, b, F, T_hat,
                                         dim, Gap_Diffraction_UInt);
        }
        else if (typenum == NPY_LONG){
            __get_complex_from_four_real(((long *)data), a, b, F, T_hat, dim,
                                         Gap_Diffraction_Long);
        }
        else if (typenum == NPY_ULONG){
            __get_complex_from_four_real(((unsigned long *)data), a, b, F,
                                         T_hat, dim, Gap_Diffraction_ULong);
        }
        else if (typenum == NPY_LONGLONG){
            __get_complex_from_four_real(((long long *)data), a, b, F, T_hat,
                                         dim, Gap_Diffraction_Long_Long);
        }
        else if (typenum == NPY_ULONG){
            __get_complex_from_four_real(((unsigned long long *)data), a, b, F,
                                         T_hat, dim,
                                         Gap_Diffraction_ULong_Long);
        }
        else {
            PyErr_Format(
                PyExc_TypeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tdiffrec.special_functions.gap_diffraction\n\n"
                "\rInvalid data type for input numpy array. Input should be\n"
                "\ra one dimensional numpy array of positive (floating point)\n"
                "\rreal numbers."
            );
            return NULL;
        }

        output = PyArray_SimpleNewFromData(1, &dim, NPY_CDOUBLE, (void *)T_hat);
        capsule = PyCapsule_New(T_hat, NULL, capsule_cleanup);
    }

    /*  This frees the variable at the Python level once it's destroyed.      */
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);

    /*  Return the results to Python.                                         */
    return Py_BuildValue("N", output);
}

static PyObject *ringlet_diffraction(PyObject *self, PyObject *args)
{
    PyObject *output, *capsule;
    PyArrayObject *rho;
    double a, b, F;
    char typenum;
    long dim;
    void *data;

    if (!PyArg_ParseTuple(args, "O!ddd", &PyArray_Type, &rho, &a, &b, &F)){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.ringlet_diffraction\n\n"
            "\rCould not parse inputs. Legal inputs are:\n"
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"
            "\r\ta:     Positive constant (Float)\n"
            "\r\tb:     Positive constant (Float) greater than a\n"
            "\r\tF      Positive constant (Float)\n\n"
            "\rNotes:\n"
            "\r\trho must be a non-empty one dimensional numpy array."
        );
        return NULL;
    }

    /*  Useful information about the data.                                */
    typenum = (char)PyArray_TYPE(rho);
    dim     = PyArray_DIMS(rho)[0];
    data    = PyArray_DATA(rho);

    /*  Check the inputs to make sure they're valid.                          */
    if (PyArray_NDIM(rho) != 1){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.ringlet_diffraction\n\n"
            "\rInput numpy array is not one-dimensional.\n"
        );
        return NULL;
    }
    else if (a >= b){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.ringlet_diffraction\n\n"
            "\rInner radius is not less than outer radius (i.e. a >= b).\n"
        );
        return NULL;
    }
    else if (a <= 0.0){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.ringlet_diffraction\n\n"
            "\rInner radius is negative. (i.e. a<0)\n"
        );
        return NULL;
    }
    else if (F < 0.0){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.ringlet_diffraction\n\n"
            "\rFresnel scale is negative (i.e. F<0).\n"
        );
        return NULL;
    }
    else if (dim == 0){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.ringlet_diffraction\n\n"
            "\rInput numpy array is empty.\n"
        );
    }

    if (typenum == NPY_FLOAT){
        complex float *T_hat;
        T_hat = (complex float *)malloc(dim*sizeof(complex float));

        __get_complex_from_four_real(
            ((float *)data), a, b, F, T_hat, dim,
            Ringlet_Diffraction_Float
        );

        output = PyArray_SimpleNewFromData(1, &dim, NPY_CFLOAT, (void *)T_hat);
        capsule = PyCapsule_New(T_hat, NULL, capsule_cleanup);
    }
    else if (typenum == NPY_LONGDOUBLE){
        complex long double *T_hat;
        T_hat = (complex long double *)malloc(dim*sizeof(complex long double));

        __get_complex_from_four_real(
            ((long double *)data), a, b, F, T_hat, dim,
            Ringlet_Diffraction_Long_Double
        );

        output = PyArray_SimpleNewFromData(1, &dim, NPY_CLONGDOUBLE,
                                           (void *)T_hat);
        capsule = PyCapsule_New(T_hat, NULL, capsule_cleanup);
    }
    else {
        complex double *T_hat;
        T_hat = (complex double *)malloc(dim*sizeof(complex double));

        if (typenum == NPY_DOUBLE){
            __get_complex_from_four_real(
                ((double *)data), a, b, F, T_hat, dim,
                Ringlet_Diffraction_Double
            );
        }
        else if (typenum == NPY_BYTE){
            __get_complex_from_four_real(((char *)data), a, b, F, T_hat, dim,
                                         Ringlet_Diffraction_Char);
        }
        else if (typenum == NPY_UBYTE){
            __get_complex_from_four_real(((unsigned char *)data), a, b, F,
                                         T_hat, dim, Ringlet_Diffraction_UChar);
        }
        else if (typenum == NPY_SHORT){
            __get_complex_from_four_real(((short *)data), a, b, F, T_hat, dim,
                                         Ringlet_Diffraction_Short);
        }
        else if (typenum == NPY_USHORT){
            __get_complex_from_four_real(((unsigned short *)data), a, b, F,
                                         T_hat, dim,
                                         Ringlet_Diffraction_UShort);
        }
        else if (typenum == NPY_INT){
            __get_complex_from_four_real(((int *)data), a, b, F, T_hat, dim,
                                         Ringlet_Diffraction_Int);
        }
        else if (typenum == NPY_UINT){
            __get_complex_from_four_real(((unsigned int *)data), a, b, F, T_hat,
                                         dim, Ringlet_Diffraction_UInt);
        }
        else if (typenum == NPY_LONG){
            __get_complex_from_four_real(((long *)data), a, b, F, T_hat, dim,
                                         Ringlet_Diffraction_Long);
        }
        else if (typenum == NPY_ULONG){
            __get_complex_from_four_real(((unsigned long *)data), a, b, F,
                                         T_hat, dim, Ringlet_Diffraction_ULong);
        }
        else if (typenum == NPY_LONGLONG){
            __get_complex_from_four_real(((long long *)data), a, b, F, T_hat,
                                         dim, Ringlet_Diffraction_Long_Long);
        }
        else if (typenum == NPY_ULONG){
            __get_complex_from_four_real(((unsigned long long *)data), a, b, F,
                                         T_hat, dim,
                                         Ringlet_Diffraction_ULong_Long);
        }
        else {
            PyErr_Format(
                PyExc_TypeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tdiffrec.special_functions.ringlet_diffraction\n\n"
                "\rInvalid data type for input numpy array. Input should be\n"
                "\ra one dimensional numpy array of positive (floating point)\n"
                "\rreal numbers."
            );
            return NULL;
        }

        output = PyArray_SimpleNewFromData(1, &dim, NPY_CDOUBLE, (void *)T_hat);
        capsule = PyCapsule_New(T_hat, NULL, capsule_cleanup);
    }

    /*  This frees the variable at the Python level once it's destroyed.      */
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);

    /*  Return the results to Python.                                         */
    return Py_BuildValue("N", output);
}

static PyObject *right_straightedge(PyObject *self, PyObject *args)
{
    PyObject *output, *capsule;
    PyArrayObject *rho;
    double a, F;
    char typenum;
    long dim;
    void *data;

    if (!PyArg_ParseTuple(args, "O!dd", &PyArray_Type, &rho, &a, &F)){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.right_straightedge\n\n"
            "\rCould not parse inputs. Legal inputs are:\n"
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"
            "\r\ta:     Positive constant (Float)\n"
            "\r\tF      Positive constant (Float)\n\n"
            "\rNotes:\n"
            "\r\trho must be a non-empty one dimensional numpy array."
        );
        return NULL;
    }

    /*  Useful information about the data.                                */
    typenum = (char)PyArray_TYPE(rho);
    dim     = PyArray_DIMS(rho)[0];
    data    = PyArray_DATA(rho);

    /*  Check the inputs to make sure they're valid.                          */
    if (PyArray_NDIM(rho) != 1){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.right_straightedge\n\n"
            "\rInput numpy array is not one-dimensional.\n"
        );
        return NULL;
    }
    else if (a <= 0.0){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.right_straightedge\n\n"
            "\rInner radius is negative. (i.e. a<0)\n"
        );
        return NULL;
    }
    else if (F < 0.0){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.right_straightedge\n\n"
            "\rFresnel scale is negative (i.e. F<0).\n"
        );
        return NULL;
    }
    else if (dim == 0){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.right_straightedge\n\n"
            "\rInput numpy array is empty.\n"
        );
    }

    if (typenum == NPY_FLOAT){
        complex float *T_hat;
        T_hat = (complex float *)malloc(dim*sizeof(complex float));

        __get_complex_from_three_real(
            ((float *)data), a, F, T_hat, dim,
            Right_Straightedge_Diffraction_Float
        );

        output = PyArray_SimpleNewFromData(1, &dim, NPY_CFLOAT, (void *)T_hat);
        capsule = PyCapsule_New(T_hat, NULL, capsule_cleanup);
    }
    else if (typenum == NPY_LONGDOUBLE){
        complex long double *T_hat;
        T_hat = (complex long double *)malloc(dim*sizeof(complex long double));

        __get_complex_from_three_real(
            ((long double *)data), a, F, T_hat, dim,
            Right_Straightedge_Diffraction_Long_Double
        );

        output = PyArray_SimpleNewFromData(1, &dim, NPY_CLONGDOUBLE,
                                           (void *)T_hat);
        capsule = PyCapsule_New(T_hat, NULL, capsule_cleanup);
    }
    else {
        complex double *T_hat;
        T_hat = (complex double *)malloc(dim*sizeof(complex double));

        if (typenum == NPY_DOUBLE){
            __get_complex_from_three_real(
                ((double *)data), a, F, T_hat, dim,
                Right_Straightedge_Diffraction_Double
            );
        }
        else if (typenum == NPY_BYTE){
            __get_complex_from_three_real(
                ((char *)data), a, F, T_hat, dim,
                Right_Straightedge_Diffraction_Char
            );
        }
        else if (typenum == NPY_UBYTE){
            __get_complex_from_three_real(
                ((unsigned char *)data), a, F, T_hat, dim,
                Right_Straightedge_Diffraction_UChar
            );
        }
        else if (typenum == NPY_SHORT){
            __get_complex_from_three_real(
                ((short *)data), a, F, T_hat, dim,
                Right_Straightedge_Diffraction_Short
            );
        }
        else if (typenum == NPY_USHORT){
            __get_complex_from_three_real(
                ((unsigned short *)data), a, F, T_hat, dim,
                Right_Straightedge_Diffraction_UShort
            );
        }
        else if (typenum == NPY_INT){
            __get_complex_from_three_real(
                ((int *)data), a, F, T_hat, dim,
                Right_Straightedge_Diffraction_Int
            );
        }
        else if (typenum == NPY_UINT){
            __get_complex_from_three_real(
                ((unsigned int *)data), a, F, T_hat, dim,
                Right_Straightedge_Diffraction_UInt);
        }
        else if (typenum == NPY_LONG){
            __get_complex_from_three_real(
                ((long *)data), a, F, T_hat, dim,
                Right_Straightedge_Diffraction_Long
            );
        }
        else if (typenum == NPY_ULONG){
            __get_complex_from_three_real(
                ((unsigned long *)data), a, F, T_hat, dim,
                Right_Straightedge_Diffraction_ULong
            );
        }
        else if (typenum == NPY_LONGLONG){
            __get_complex_from_three_real(
                ((long long *)data), a, F, T_hat, dim,
                Right_Straightedge_Diffraction_Long_Long);
        }
        else if (typenum == NPY_ULONG){
            __get_complex_from_three_real(
                ((unsigned long long *)data), a, F, T_hat, dim,
                Right_Straightedge_Diffraction_ULong_Long);
        }
        else {
            PyErr_Format(
                PyExc_TypeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tdiffrec.special_functions.right_straightedge\n\n"
                "\rInvalid data type for input numpy array. Input should be\n"
                "\ra one dimensional numpy array of positive (floating point)\n"
                "\rreal numbers."
            );
            return NULL;
        }

        output = PyArray_SimpleNewFromData(1, &dim, NPY_CDOUBLE, (void *)T_hat);
        capsule = PyCapsule_New(T_hat, NULL, capsule_cleanup);
    }

    /*  This frees the variable at the Python level once it's destroyed.      */
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);

    /*  Return the results to Python.                                         */
    return Py_BuildValue("N", output);
}

static PyObject *left_straightedge(PyObject *self, PyObject *args)
{
    PyObject *output, *capsule;
    PyArrayObject *rho;
    double a, F;
    char typenum;
    long dim;
    void *data;

    if (!PyArg_ParseTuple(args, "O!dd", &PyArray_Type, &rho, &a, &F)){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.left_straightedge\n\n"
            "\rCould not parse inputs. Legal inputs are:\n"
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"
            "\r\ta:     Positive constant (Float)\n"
            "\r\tF      Positive constant (Float)\n\n"
            "\rNotes:\n"
            "\r\trho must be a non-empty one dimensional numpy array."
        );
        return NULL;
    }

    /*  Useful information about the data.                                */
    typenum = (char)PyArray_TYPE(rho);
    dim     = PyArray_DIMS(rho)[0];
    data    = PyArray_DATA(rho);

    /*  Check the inputs to make sure they're valid.                          */
    if (PyArray_NDIM(rho) != 1){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.left_straightedge\n\n"
            "\rInput numpy array is not one-dimensional.\n"
        );
        return NULL;
    }
    else if (a <= 0.0){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.left_straightedge\n\n"
            "\rInner radius is negative. (i.e. a<0)\n"
        );
        return NULL;
    }
    else if (F < 0.0){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.left_straightedge\n\n"
            "\rFresnel scale is negative (i.e. F<0).\n"
        );
        return NULL;
    }
    else if (dim == 0){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.left_straightedge\n\n"
            "\rInput numpy array is empty.\n"
        );
    }

    if (typenum == NPY_FLOAT){
        complex float *T_hat;
        T_hat = (complex float *)malloc(dim*sizeof(complex float));

        __get_complex_from_three_real(
            ((float *)data), a, F, T_hat, dim,
            Left_Straightedge_Diffraction_Float
        );

        output = PyArray_SimpleNewFromData(1, &dim, NPY_CFLOAT, (void *)T_hat);
        capsule = PyCapsule_New(T_hat, NULL, capsule_cleanup);
    }
    else if (typenum == NPY_LONGDOUBLE){
        complex long double *T_hat;
        T_hat = (complex long double *)malloc(dim*sizeof(complex long double));

        __get_complex_from_three_real(
            ((long double *)data), a, F, T_hat, dim,
            Left_Straightedge_Diffraction_Long_Double
        );

        output = PyArray_SimpleNewFromData(1, &dim, NPY_CLONGDOUBLE,
                                           (void *)T_hat);
        capsule = PyCapsule_New(T_hat, NULL, capsule_cleanup);
    }
    else {
        complex double *T_hat;
        T_hat = (complex double *)malloc(dim*sizeof(complex double));

        if (typenum == NPY_DOUBLE){
            __get_complex_from_three_real(
                ((double *)data), a, F, T_hat, dim,
                Left_Straightedge_Diffraction_Double
            );
        }
        else if (typenum == NPY_BYTE){
            __get_complex_from_three_real(
                ((char *)data), a, F, T_hat, dim,
                Left_Straightedge_Diffraction_Char
            );
        }
        else if (typenum == NPY_UBYTE){
            __get_complex_from_three_real(
                ((unsigned char *)data), a, F, T_hat, dim,
                Left_Straightedge_Diffraction_UChar
            );
        }
        else if (typenum == NPY_SHORT){
            __get_complex_from_three_real(
                ((short *)data), a, F, T_hat, dim,
                Left_Straightedge_Diffraction_Short
            );
        }
        else if (typenum == NPY_USHORT){
            __get_complex_from_three_real(
                ((unsigned short *)data), a, F, T_hat, dim,
                Left_Straightedge_Diffraction_UShort
            );
        }
        else if (typenum == NPY_INT){
            __get_complex_from_three_real(
                ((int *)data), a, F, T_hat, dim,
                Left_Straightedge_Diffraction_Int
            );
        }
        else if (typenum == NPY_UINT){
            __get_complex_from_three_real(
                ((unsigned int *)data), a, F, T_hat, dim,
                Left_Straightedge_Diffraction_UInt);
        }
        else if (typenum == NPY_LONG){
            __get_complex_from_three_real(
                ((long *)data), a, F, T_hat, dim,
                Left_Straightedge_Diffraction_Long
            );
        }
        else if (typenum == NPY_ULONG){
            __get_complex_from_three_real(
                ((unsigned long *)data), a, F, T_hat, dim,
                Left_Straightedge_Diffraction_ULong
            );
        }
        else if (typenum == NPY_LONGLONG){
            __get_complex_from_three_real(
                ((long long *)data), a, F, T_hat, dim,
                Left_Straightedge_Diffraction_Long_Long);
        }
        else if (typenum == NPY_ULONG){
            __get_complex_from_three_real(
                ((unsigned long long *)data), a, F, T_hat, dim,
                Left_Straightedge_Diffraction_ULong_Long);
        }
        else {
            PyErr_Format(
                PyExc_TypeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tdiffrec.special_functions.left_straightedge\n\n"
                "\rInvalid data type for input numpy array. Input should be\n"
                "\ra one dimensional numpy array of positive (floating point)\n"
                "\rreal numbers."
            );
            return NULL;
        }

        output = PyArray_SimpleNewFromData(1, &dim, NPY_CDOUBLE, (void *)T_hat);
        capsule = PyCapsule_New(T_hat, NULL, capsule_cleanup);
    }

    /*  This frees the variable at the Python level once it's destroyed.      */
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);

    /*  Return the results to Python.                                         */
    return Py_BuildValue("N", output);
}

static PyObject *square_wave_diffraction(PyObject *self, PyObject *args)
{
    PyObject *output, *capsule;
    PyArrayObject *x_arr;
    double W, F;
    long N;
    char typenum;
    long dim;
    void *data;
    complex double *T_hat;

    if (!PyArg_ParseTuple(args, "O!ddK", &PyArray_Type, &x_arr, &W, &F, &N)){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.square_wave_diffraction\n\n"
            "\rCould not parse inputs. Legal inputs are:\n"
            "\r\tx:     Numpy Array of positive real numbers (Floats)\n"
            "\r\tW:     Positive constant (Float)\n"
            "\r\tF:     Positive constant (Float)\n"
            "\r\tN:     Positive Integer (Int)\n\n"
            "\rNotes:\n"
            "\r\trho must be a non-empty one dimensional numpy array."
        );
        return NULL;
    }

    /*  Useful information about the data.                                */
    typenum = (char)PyArray_TYPE(x_arr);
    dim     = PyArray_DIMS(x_arr)[0];
    data    = PyArray_DATA(x_arr);

    /*  Check the inputs to make sure they're valid.                          */
    if (PyArray_NDIM(x_arr) != 1){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.square_wave_diffraction\n\n"
            "\rInput numpy array is not one-dimensional.\n"
        );
        return NULL;
    }
    else if (W <= 0.0){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.square_wave_diffraction\n\n"
            "\rWidth of wave is non-positive. (i.e. W<=0)\n"
        );
        return NULL;
    }
    else if (F <= 0.0){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.square_wave_diffraction\n\n"
            "\rFresnel scale is non-positive (i.e. F<=0).\n"
        );
        return NULL;
    }
    else if (dim == 0){
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.square_wave_diffraction\n\n"
            "\rInput numpy array is empty.\n"
        );
    }

    T_hat = (complex double *)malloc(dim*sizeof(complex double));
    if (typenum == NPY_DOUBLE){
            __get_complex_from_four_real(
                ((double *)data), W, F, N, T_hat, dim,
                Square_Wave_Diffraction_Double
            );
        }
    else {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.special_functions.square_wave_diffraction\n\n"
            "\rInvalid data type for input numpy array. Input should be\n"
            "\ra one dimensional numpy array of positive (floating point)\n"
            "\rreal numbers."
        );
        return NULL;
    }

    output = PyArray_SimpleNewFromData(1, &dim, NPY_CDOUBLE, (void *)T_hat);
    capsule = PyCapsule_New(T_hat, NULL, capsule_cleanup);

    /*  This frees the variable at the Python level once it's destroyed.      */
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);

    /*  Return the results to Python.                                         */
    return Py_BuildValue("N", output);
}


static PyMethodDef _special_functions_methods[] =
{
    {
        "compute_norm_eq",
        compute_norm_eq,
        METH_VARARGS,
        "Compute the normalized equivalent width of an array."
    },
    {
        "fresnel_transform",
        fresnel_transform,
        METH_VARARGS,
        "Compute the Fresnel Transform."
    },
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
        "_special_functions.gap_diffraction\n\r\t"
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
        ">>> import _special_functions\n\r\t\t"
        ">>> x = numpy.arange(0,100,0.01)\n\r\t\t"
        ">>> a = 45\n\r\t\t"
        ">>> b = 55\n\r\t\t"
        ">>> F = 0.05\n\r\t\t"
        ">>> y = _special_functions.gap_diffraction(x, 45, 55, 0.05)"
    },
    {
        "ringlet_diffraction",
        ringlet_diffraction,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "_special_functions.gap_diffraction\n\r\t"
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
        ">>> import _special_functions\n\r\t\t"
        ">>> x = numpy.arange(0,100,0.01)\n\r\t\t"
        ">>> a = 45\n\r\t\t"
        ">>> b = 55\n\r\t\t"
        ">>> F = 0.05\n\r\t\t"
        ">>> y = _special_functions.ringlet_diffraction(x, a, b, F)"
    },
    {
        "right_straightedge",
        right_straightedge,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "_special_functions.right_straightedge\n\r\t"
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
        ">>> import _special_functions\n\r\t\t"
        ">>> x = numpy.arange(0,100,0.01)\n\r\t\t"
        ">>> a = 45\n\r\t\t"
        ">>> F = 0.05\n\r\t\t"
        ">>> y = _special_functions.right_straightedge(x, a, F)"
    },
    {
        "left_straightedge",
        left_straightedge,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "_special_functions.left_straightedge\n\r\t"
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
        ">>> import _special_functions\n\r\t\t"
        ">>> x = numpy.arange(0,100,0.01)\n\r\t\t"
        ">>> a = 45\n\r\t\t"
        ">>> F = 0.05\n\r\t\t"
        ">>> y = _special_functions.left_straightedge(x, a, F)"
    },
    {
        "max",
        max,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "_special_functions.max\n\r\t"
        "Purpose\n\r\t\t"
        "Compute the maximum of a numpy array.\n\r\t"
        "Arguments:\n\r\t\t"
        "arr (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers.\n\r\t"
        "Outputs:\n\r\t\t"
        "max (int or float):\n\r\t\t\t"
        "The maximum value of the input array.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import _special_functions\n\r\t\t"
        ">>> x = numpy.random.rand(100)\n\r\t\t"
        ">>> y = _special_functions.max(x)"
    },
    {
        "min",
        min,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "_special_functions.min\n\r\t"
        "Purpose\n\r\t\t"
        "Compute the minimum of a numpy array.\n\r\t"
        "Arguments:\n\r\t\t"
        "arr (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers.\n\r\t"
        "Outputs:\n\r\t\t"
        "min (int or float):\n\r\t\t\t"
        "The minimum value of the input array.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import _special_functions\n\r\t\t"
        ">>> x = numpy.random.rand(100)\n\r\t\t"
        ">>> y = _special_functions.min(x)"
    },
    {
        "where_greater",
        where_greater,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "_special_functions.where_greater\n\r\t"
        "Purpose\n\r\t\t"
        "Given a real-valued numpy array arr, and a real number\n\r\t\t"
        "threshold, compute the indices n such that arr[n] > threshold\n\r\t"
        "Arguments:\n\r\t\t"
        "arr (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers.\n\r\t\t"
        "threshold (int or float):\n\r\t\t\t"
        "The threshold value for comparing arr with."
        "Outputs:\n\r\t\t"
        "where_arr (numpy.ndarray):\n\r\t\t\t"
        "The array of indices such that arr[n] > threshold.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import _special_functions\n\r\t\t"
        ">>> x = numpy.arange(-5, 5, 0.01)\n\r\t\t"
        ">>> y = _special_functions.where_greater(x, 1.0)"
    },
    {
        "where_lesser",
        where_lesser,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "_special_functions.where_lesser\n\r\t"
        "Purpose\n\r\t\t"
        "Given a real-valued numpy array arr, and a real number\n\r\t\t"
        "threshold, compute the indices n such that arr[n] < threshold\n\r\t"
        "Arguments:\n\r\t\t"
        "arr (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of real numbers.\n\r\t\t"
        "threshold (int or float):\n\r\t\t\t"
        "The threshold value for comparing arr with."
        "Outputs:\n\r\t\t"
        "where_arr (numpy.ndarray):\n\r\t\t\t"
        "The array of indices such that arr[n] < threshold.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import _special_functions\n\r\t\t"
        ">>> x = numpy.arange(-5, 5, 0.01)\n\r\t\t"
        ">>> y = _special_functions.where_lesser(x, 1.0)"
    },
    {NULL, NULL, 0, NULL}
};
/*-------------------------DEFINE UNIVERSAL FUNCTIONS-------------------------*/

/*------------------------------BESSEL J0-------------------------------------*/
static void char_besselJ0(char **args, npy_intp *dimensions,
                          npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    char *x   = (char *)args[0];
    double *y = (double *)args[1];

    Get_Char_to_Double_Array(x, y, n_elements, BesselJ0_Char);
}

static void uchar_besselJ0(char **args, npy_intp *dimensions,
                           npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    unsigned char *x = (unsigned char *)args[0];
    double *y        = (double *)args[1];

    Get_UChar_to_Double_Array(x, y, n_elements, BesselJ0_UChar);
}

static void short_besselJ0(char **args, npy_intp *dimensions,
                           npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    short *x  = (short *)args[0];
    double *y = (double *)args[1];

    Get_Short_to_Double_Array(x, y, n_elements, BesselJ0_Short);
}

static void ushort_besselJ0(char **args, npy_intp *dimensions,
                            npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    unsigned short *x = (unsigned short *)args[0];
    double *y         = (double *)args[1];

    Get_UShort_to_Double_Array(x, y, n_elements, BesselJ0_UShort);
}

static void int_besselJ0(char **args, npy_intp *dimensions,
                         npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    int *x    = (int *)args[0];
    double *y = (double *)args[1];

    Get_Int_to_Double_Array(x, y, n_elements, BesselJ0_Int);
}

static void uint_besselJ0(char **args, npy_intp *dimensions,
                          npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    unsigned int *x = (unsigned int *)args[0];
    double *y       = (double *)args[1];

    Get_UInt_to_Double_Array(x, y, n_elements, BesselJ0_UInt);
}

static void long_besselJ0(char **args, npy_intp *dimensions,
                          npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    long *x   = (long *)args[0];
    double *y = (double *)args[1];

    Get_Long_to_Double_Array(x, y, n_elements, BesselJ0_Long);
}

static void ulong_besselJ0(char **args, npy_intp *dimensions,
                           npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    unsigned long *x = (unsigned long *)args[0];
    double *y        = (double *)args[1];

    Get_ULong_to_Double_Array(x, y, n_elements, BesselJ0_ULong);
}

static void long_long_besselJ0(char **args, npy_intp *dimensions,
                               npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    long long *x = (long long *)args[0];
    double *y    = (double *)args[1];

    Get_Long_Long_to_Double_Array(x, y, n_elements, BesselJ0_Long_Long);
}

static void ulong_long_besselJ0(char **args, npy_intp *dimensions,
                                npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    unsigned long long *x = (unsigned long long *)args[0];
    double *y             = (double *)args[1];

    Get_ULong_Long_to_Double_Array(x, y, n_elements, BesselJ0_ULong_Long);
}

static void float_besselJ0(char **args, npy_intp *dimensions,
                           npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    float *x = (float *)args[0];
    float *y = (float *)args[1];

    Get_Float_Array(x, y, n_elements, BesselJ0_Float);
}

static void double_besselJ0(char **args, npy_intp *dimensions,
                            npy_intp *steps, void *data)
{

    npy_intp n_elements = dimensions[0];

    double *x = (double *)args[0];
    double *y = (double *)args[1];

    Get_Double_Array(x, y, n_elements, BesselJ0_Double);
}

static void long_double_besselJ0(char **args, npy_intp *dimensions,
                                 npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    long double *x = (long double *)args[0];
    long double *y = (long double *)args[1];

    Get_Long_Double_Array(x, y, n_elements, BesselJ0_Long_Double);
}


PyUFuncGenericFunction besselJ0_funcs[13] = {
    &char_besselJ0,
    &uchar_besselJ0,
    &short_besselJ0,
    &ushort_besselJ0,
    &int_besselJ0,
    &uint_besselJ0,
    &long_besselJ0,
    &ulong_besselJ0,
    &long_long_besselJ0,
    &ulong_long_besselJ0,
    &float_besselJ0,
    &double_besselJ0,
    &long_double_besselJ0
};

/*------------------------------BESSEL I0-------------------------------------*/
static void char_besselI0(char **args, npy_intp *dimensions,
                          npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    char *x   = (char *)args[0];
    double *y = (double *)args[1];

    Get_Char_to_Double_Array(x, y, n_elements, BesselI0_Char);
}

static void uchar_besselI0(char **args, npy_intp *dimensions,
                           npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    unsigned char *x = (unsigned char *)args[0];
    double *y        = (double *)args[1];

    Get_UChar_to_Double_Array(x, y, n_elements, BesselI0_UChar);
}

static void short_besselI0(char **args, npy_intp *dimensions,
                           npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    short *x  = (short *)args[0];
    double *y = (double *)args[1];

    Get_Short_to_Double_Array(x, y, n_elements, BesselI0_Short);
}

static void ushort_besselI0(char **args, npy_intp *dimensions,
                            npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    unsigned short *x = (unsigned short *)args[0];
    double *y         = (double *)args[1];

    Get_UShort_to_Double_Array(x, y, n_elements, BesselI0_UShort);
}

static void int_besselI0(char **args, npy_intp *dimensions,
                         npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    int *x    = (int *)args[0];
    double *y = (double *)args[1];

    Get_Int_to_Double_Array(x, y, n_elements, BesselI0_Int);
}

static void uint_besselI0(char **args, npy_intp *dimensions,
                          npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    unsigned int *x = (unsigned int *)args[0];
    double *y       = (double *)args[1];

    Get_UInt_to_Double_Array(x, y, n_elements, BesselI0_UInt);
}

static void long_besselI0(char **args, npy_intp *dimensions,
                          npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    long *x   = (long *)args[0];
    double *y = (double *)args[1];

    Get_Long_to_Double_Array(x, y, n_elements, BesselI0_Long);
}

static void ulong_besselI0(char **args, npy_intp *dimensions,
                           npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    unsigned long *x = (unsigned long *)args[0];
    double *y        = (double *)args[1];

    Get_ULong_to_Double_Array(x, y, n_elements, BesselI0_ULong);
}

static void long_long_besselI0(char **args, npy_intp *dimensions,
                               npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    long long *x = (long long *)args[0];
    double *y    = (double *)args[1];

    Get_Long_Long_to_Double_Array(x, y, n_elements, BesselI0_Long_Long);
}

static void ulong_long_besselI0(char **args, npy_intp *dimensions,
                                npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    unsigned long long *x = (unsigned long long *)args[0];
    double *y             = (double *)args[1];

    Get_ULong_Long_to_Double_Array(x, y, n_elements, BesselI0_ULong_Long);
}

static void float_besselI0(char **args, npy_intp *dimensions,
                           npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    float *x = (float *)args[0];
    float *y = (float *)args[1];

    Get_Float_Array(x, y, n_elements, BesselI0_Float);
}

static void double_besselI0(char **args, npy_intp *dimensions,
                            npy_intp *steps, void *data)
{

    npy_intp n_elements = dimensions[0];

    double *x = (double *)args[0];
    double *y = (double *)args[1];

    Get_Double_Array(x, y, n_elements, BesselI0_Double);
}

static void long_double_besselI0(char **args, npy_intp *dimensions,
                                 npy_intp *steps, void *data)
{
    npy_intp n_elements = dimensions[0];

    long double *x = (long double *)args[0];
    long double *y = (long double *)args[1];

    Get_Long_Double_Array(x, y, n_elements, BesselI0_Long_Double);
}


PyUFuncGenericFunction besselI0_funcs[13] = {
    &char_besselI0,
    &uchar_besselI0,
    &short_besselI0,
    &ushort_besselI0,
    &int_besselI0,
    &uint_besselI0,
    &long_besselI0,
    &ulong_besselI0,
    &long_long_besselI0,
    &ulong_long_besselI0,
    &float_besselI0,
    &double_besselI0,
    &long_double_besselI0
};

/************Double Slit Diffraction Using Fraunhofer Approximation************/

/******************************************************************************
 *  Function:                                                                 *
 *      float_double_slit_diffraction                                         *
 *  Purpose:                                                                  *
 *      Compute the diffraction pattern from a plane wave incident on a       *
 *      single slit using the Fraunhofer approximation.                       *
 *  Arguments:                                                                *
 *      args (char **):                                                       *
 *          Input and output arguments passed from python.                    *
 *      dimensions (npy_intp *):                                              *
 *          Dimensions of the arguments found in the args pointer.            *
 *      steps (npy_intp):                                                     *
 *          The number of strides in memory from the nth point to the (n+1)th *
 *          point for the arguments found in the args pointer.                *
 *      data (void *):                                                        *
 *          Data pointer.                                                     *
 *  Notes:                                                                    *
 *      1.) This is a wrapper for Double_Slit_Fraunhofer_Diffraction_Float,   *
 *          which is defined in __fraunhofer_diffraction.h. This allows       *
 *          Python to use that function, and allows for numpy arrays to be    *
 *          passed in. Relies on the Numpy UFUNC API and the C-Python API.    *
 *                                                                            *
 *      2.) This function relies on the C99 standard, or higher.              *
 *                                                                            *
 *      3.) There are no error checks in this code. This is handled at the    *
 *          Python level, see special_functions.py.                           *
 ******************************************************************************/
static void float_double_slit_diffraction(char **args, npy_intp *dimensions,
                                          npy_intp* steps, void* data)
{
    /* Declare i for indexing, n is the number of elements in the array.      */
    npy_intp i;
    npy_intp n = dimensions[0];

    /* Extract input data and convert to appropriate types.                   */
    float *x  =  (float *)args[0];
    float  z  = *(float *)args[1];
    float  a  = *(float *)args[2];
    float  d  = *(float *)args[3];
    float  l  = *(float *)args[4];

    /* The output is a pointer to a complex float.                            */
    float *out = (float *)args[5];

    /* Loop over the square well function found in __fresnel_diffraction.h    */
    for (i = 0; i < n; i++) {
        out[i] = Double_Slit_Fraunhofer_Diffraction_Float(x[i], z, a, d, l);
    }
}

static void double_double_slit_diffraction(char **args, npy_intp *dimensions,
                                           npy_intp* steps, void* data)
{
    /* Declare i for indexing, n is the number of elements in the array.      */
    npy_intp i;
    npy_intp n = dimensions[0];

    /* Extract input data and convert to appropriate types.                   */
    double *x  =  (double *)args[0];
    double  z  = *(double *)args[1];
    double  a  = *(double *)args[2];
    double  d  = *(double *)args[3];
    double  l  = *(double *)args[4];

    /* The output is a pointer to a complex float.                            */
    double *out = (double *)args[5];

    /* Loop over the square well function found in __fresnel_diffraction.h    */
    for (i = 0; i < n; i++) {
        out[i] = Double_Slit_Fraunhofer_Diffraction_Double(x[i], z, a, d, l);
    }
}

static void long_double_double_slit_diffraction(char **args,
                                                npy_intp *dimensions,
                                                npy_intp* steps, void* data)
{
    /* Declare i for indexing, n is the number of elements in the array.      */
    npy_intp i;
    npy_intp n = dimensions[0];

    /* Extract input data and convert to appropriate types.                   */
    long double *x  =  (long double *)args[0];
    long double  z  = *(long double *)args[1];
    long double  a  = *(long double *)args[2];
    long double  d  = *(long double *)args[3];
    long double  l  = *(long double *)args[4];

    /* The output is a pointer to a complex float.                            */
    long double *out = (long double *)args[5];

    /* Loop over the square well function found in __fresnel_diffraction.h    */
    for (i = 0; i < n; i++) {
        out[i] = Double_Slit_Fraunhofer_Diffraction_Long_Double(x[i], z,
                                                                a, d, l);
    }
}

PyUFuncGenericFunction double_slit_funcs[3] = {
    &float_double_slit_diffraction,
    &double_double_slit_diffraction,
    &long_double_double_slit_diffraction
};

PyUFuncGenericFunction frequency_to_wavelength_funcs[3] = {
    &float_frequency_to_wavelength,
    &double_frequency_to_wavelength,
    &long_double_frequency_to_wavelength
};

PyUFuncGenericFunction fresnel_cos_funcs[3] = {
    &float_fresnelcos,
    &double_fresnelcos,
    &long_double_fresnelcos

};

PyUFuncGenericFunction fresnel_scale_funcs[3] = {
    &float_fresnel_scale,
    &double_fresnel_scale,
    &long_double_fresnel_scale
};

PyUFuncGenericFunction fresnel_sin_funcs[3] = {
    &float_fresnelsin,
    &double_fresnelsin,
    &long_double_fresnelsin
};

PyUFuncGenericFunction lambertw_funcs[3] = {
    &float_lambertw,
    &double_lambertw,
    &long_double_lambertw
};

PyUFuncGenericFunction psi_funcs[3] = {
    &float_fresnel_psi,
    &double_fresnel_psi,
    &long_double_fresnel_psi
};

PyUFuncGenericFunction dpsi_funcs[3] = {
    &float_fresnel_dpsi_dphi,
    &double_fresnel_dpsi_dphi,
    &long_double_fresnel_dpsi_dphi
};

PyUFuncGenericFunction d2psi_funcs[3] = {
    &float_fresnel_d2psi_dphi2,
    &double_fresnel_d2psi_dphi2,
    &long_double_fresnel_d2psi_dphi2
};

PyUFuncGenericFunction dpsi_ellipse_funcs[3] = {
    &float_fresnel_dpsi_dphi_ellipse,
    &double_fresnel_dpsi_dphi_ellipse,
    &long_double_fresnel_dpsi_dphi_ellipse
};

PyUFuncGenericFunction res_inv_funcs[3] = {
    &float_resolution_inverse,
    &double_resolution_inverse,
    &long_double_resolution_inverse
};

PyUFuncGenericFunction sinc_funcs[3] = {
    &float_sinc,
    &double_sinc,
    &long_double_sinc
};

PyUFuncGenericFunction single_slit_funcs[3] = {
    &float_single_slit_diffraction,
    &double_single_slit_diffraction,
    &long_double_single_slit_diffraction
};

PyUFuncGenericFunction wavelength_to_wavenumber_funcs[3] = {
    &float_wavelength_to_wavenumber,
    &double_wavelength_to_wavenumber,
    &long_double_wavelength_to_wavenumber
};

/*  Input and return types for double input and out.                          */
static void *PyuFunc_None_3[3] = {NULL, NULL, NULL};
static void *PyuFunc_None_13[13] = {
    NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL
};

static char one_real_in_one_real_out[6] = {NPY_FLOAT, NPY_FLOAT,
                                           NPY_DOUBLE, NPY_DOUBLE,
                                           NPY_LONGDOUBLE, NPY_LONGDOUBLE};

static char one_real_in_one_float_out[26] = {
    NPY_BYTE, NPY_DOUBLE,
    NPY_UBYTE, NPY_DOUBLE,
    NPY_SHORT, NPY_DOUBLE,
    NPY_USHORT, NPY_DOUBLE,
    NPY_INT, NPY_DOUBLE,
    NPY_UINT, NPY_DOUBLE,
    NPY_LONG, NPY_DOUBLE,
    NPY_ULONG, NPY_DOUBLE,
    NPY_LONGLONG, NPY_DOUBLE,
    NPY_ULONGLONG, NPY_DOUBLE,
    NPY_FLOAT, NPY_FLOAT,
    NPY_DOUBLE, NPY_DOUBLE,
    NPY_LONGDOUBLE, NPY_LONGDOUBLE
};

static char three_real_in_one_real_out[12] = {
    NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
    NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE
};

static char four_real_in_one_real_out[15] = {
    NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
    NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE,
    NPY_LONGDOUBLE
};

static char five_real_in_one_real_out[18] = {
    NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
    NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE,
    NPY_LONGDOUBLE, NPY_LONGDOUBLE
};

static char seven_real_in_one_real_out[24] = {
    NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
    NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
    NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE,
    NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE
};

static char nine_real_in_one_real_out[30] = {
    NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
    NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
    NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE,
    NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE,
    NPY_LONGDOUBLE, NPY_LONGDOUBLE
};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_special_functions",
    NULL,
    -1,
    _special_functions_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit__special_functions(void)
{
    PyObject *besselJ0;
    PyObject *besselI0;
    PyObject *double_slit_diffraction;
    PyObject *frequency_to_wavelength;
    PyObject *fresnel_cos;
    PyObject *fresnel_psi;
    PyObject *fresnel_dpsi_dphi;
    PyObject *fresnel_d2psi_dphi2;
    PyObject *fresnel_dpsi_dphi_ellipse;
    PyObject *fresnel_scale;
    PyObject *fresnel_sin;
    PyObject *lambertw;
    PyObject *resolution_inverse;
    PyObject *single_slit_diffraction;
    PyObject *sinc;
    PyObject *wavelength_to_wavenumber;
    PyObject *m, *d;

    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    import_array();
    import_umath();

    besselJ0 = PyUFunc_FromFuncAndData(
        besselJ0_funcs, PyuFunc_None_13, one_real_in_one_float_out, 13, 1, 1,
        PyUFunc_None, "besselJ0_diffraction", "besselJ0_docstring", 0
    );

    besselI0 = PyUFunc_FromFuncAndData(
        besselI0_funcs, PyuFunc_None_13, one_real_in_one_float_out, 13, 1, 1,
        PyUFunc_None, "besselI0_diffraction", "besselI0_docstring", 0
    );

    double_slit_diffraction = PyUFunc_FromFuncAndData(
        double_slit_funcs, PyuFunc_None_3, five_real_in_one_real_out,
        3, 5, 1, PyUFunc_None, "double_slit_diffraction", 
        "double_slit_diffraction_docstring", 0
    );

    frequency_to_wavelength = PyUFunc_FromFuncAndData(
        frequency_to_wavelength_funcs, PyuFunc_None_3, one_real_in_one_real_out,
        3, 1, 1, PyUFunc_None, "frequency_to_wavelength",
        "frequency_to_wavelength_docstring", 0
    );

    fresnel_dpsi_dphi = PyUFunc_FromFuncAndData(
        dpsi_funcs, PyuFunc_None_3, seven_real_in_one_real_out, 3, 7, 1,
        PyUFunc_None, "fresnel_dpsi_dphi",  "fresnel_dpsi_dphi_docstring", 0
    );

    fresnel_d2psi_dphi2 = PyUFunc_FromFuncAndData(
        d2psi_funcs, PyuFunc_None_3, seven_real_in_one_real_out, 3, 7, 1,
        PyUFunc_None, "fresnel_d2psi_dphi2",  "fresnel_d2psi_dphi2_docstring", 0
    );

    fresnel_dpsi_dphi_ellipse = PyUFunc_FromFuncAndData(
        dpsi_ellipse_funcs, PyuFunc_None_3, nine_real_in_one_real_out, 3, 9, 1,
        PyUFunc_None, "fresnel_dpsi_dphi_ellipse", 
        "fresnel_dpsi_dphi_ellipse_docstring", 0
    );

    fresnel_cos = PyUFunc_FromFuncAndData(
        fresnel_cos_funcs, PyuFunc_None_3, one_real_in_one_real_out, 3, 1, 1,
        PyUFunc_None, "fresnel_cos", "fresnel_cos_docstring", 0
    );

    fresnel_psi = PyUFunc_FromFuncAndData(
        psi_funcs, PyuFunc_None_3, seven_real_in_one_real_out, 3, 7, 1,
        PyUFunc_None, "fresnel_psi",  "fresnel_psi_docstring", 0
    );

    fresnel_scale = PyUFunc_FromFuncAndData(
        fresnel_scale_funcs, PyuFunc_None_3, four_real_in_one_real_out, 3, 4, 1,
        PyUFunc_None, "fresnel_scale", "fresnel_scale_docstring", 0
    );

    fresnel_sin = PyUFunc_FromFuncAndData(
        fresnel_sin_funcs, PyuFunc_None_3, one_real_in_one_real_out, 3, 1, 1,
        PyUFunc_None, "fresnel_sin", "fresnel_sin_docstring", 0
    );

    lambertw = PyUFunc_FromFuncAndData(
        lambertw_funcs, PyuFunc_None_3, one_real_in_one_real_out, 3, 1, 1,
        PyUFunc_None, "lambertw", "lambertw_docstring", 0
    );

    resolution_inverse = PyUFunc_FromFuncAndData(
        res_inv_funcs, PyuFunc_None_3, one_real_in_one_real_out,
        3, 1, 1, PyUFunc_None, "resolution_inverse", 
        "resolution_inverse_docstring", 0
    );

    sinc = PyUFunc_FromFuncAndData(
        sinc_funcs, PyuFunc_None_3, one_real_in_one_real_out,
        3, 1, 1, PyUFunc_None, "sinc",  "sinc_docstring", 0
    );

    single_slit_diffraction = PyUFunc_FromFuncAndData(
        single_slit_funcs, PyuFunc_None_3, three_real_in_one_real_out,
        3, 3, 1, PyUFunc_None, "single_slit_diffraction", 
        "single_slit_diffraction_docstring", 0
    );

    wavelength_to_wavenumber = PyUFunc_FromFuncAndData(
        wavelength_to_wavenumber_funcs, PyuFunc_None_3,
        one_real_in_one_real_out, 3, 1, 1, PyUFunc_None,
        "wavelength_to_wavenumber", "wavelength_to_wavenumber_docstring", 0
    );

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "besselJ0", besselJ0);
    PyDict_SetItemString(d, "besselI0", besselI0);
    PyDict_SetItemString(d, "double_slit_diffraction", double_slit_diffraction);
    PyDict_SetItemString(d, "frequency_to_wavelength", frequency_to_wavelength);
    PyDict_SetItemString(d, "fresnel_cos", fresnel_cos);
    PyDict_SetItemString(d, "fresnel_psi", fresnel_psi);
    PyDict_SetItemString(d, "fresnel_dpsi_dphi", fresnel_dpsi_dphi);
    PyDict_SetItemString(d, "fresnel_d2psi_dphi2", fresnel_d2psi_dphi2);
    PyDict_SetItemString(d, "fresnel_dpsi_dphi_ellipse",
                         fresnel_dpsi_dphi_ellipse);
    PyDict_SetItemString(d, "fresnel_scale", fresnel_scale);
    PyDict_SetItemString(d, "fresnel_sin", fresnel_sin);
    PyDict_SetItemString(d, "lambertw", lambertw);
    PyDict_SetItemString(d, "resolution_inverse", resolution_inverse);
    PyDict_SetItemString(d, "sinc", sinc);
    PyDict_SetItemString(d, "single_slit_diffraction", single_slit_diffraction);
    PyDict_SetItemString(d, "wavelength_to_wavenumber",
                         wavelength_to_wavenumber);

    Py_DECREF(double_slit_diffraction);
    Py_DECREF(besselJ0);
    Py_DECREF(besselI0);
    Py_DECREF(frequency_to_wavelength);
    Py_DECREF(fresnel_cos);
    Py_DECREF(fresnel_psi);
    Py_DECREF(fresnel_dpsi_dphi);
    Py_DECREF(fresnel_d2psi_dphi2);
    Py_DECREF(fresnel_dpsi_dphi_ellipse);
    Py_DECREF(fresnel_scale);
    Py_DECREF(fresnel_sin);
    Py_DECREF(lambertw);
    Py_DECREF(resolution_inverse);
    Py_DECREF(sinc);
    Py_DECREF(single_slit_diffraction);
    Py_DECREF(wavelength_to_wavenumber);

    return m;
}
#else
PyMODINIT_FUNC init__funcs(void)
{
    PyObject *besselJ0;
    PyObject *besselI0;
    PyObject *double_slit_diffraction;
    PyObject *frequency_to_wavelength;
    PyObject *fresnel_cos;
    PyObject *fresnel_psi;
    PyObject *fresnel_dpsi_dphi;
    PyObject *fresnel_d2psi_dphi2;
    PyObject *fresnel_dpsi_dphi_ellipse;
    PyObject *fresnel_scale;
    PyObject *fresnel_sin;
    PyObject *lambertw;
    PyObject *resolution_inverse;
    PyObject *single_slit_diffraction;
    PyObject *sinc;
    PyObject *wavelength_to_wavenumber;
    PyObject *m, *d;

    m = Py_InitModule("__funcs", _special_functions_methods);
    if (m == NULL) {
        return;
    }

    import_array();
    import_umath();

    besselJ0 = PyUFunc_FromFuncAndData(
        besselJ0_funcs, PyuFunc_None_13, one_real_in_one_float_out, 13, 1, 1,
        PyUFunc_None, "besselJ0_diffraction", "besselJ0_docstring", 0
    );

    besselI0 = PyUFunc_FromFuncAndData(
        besselI0_funcs, PyuFunc_None_13, one_real_in_one_float_out, 13, 1, 1,
        PyUFunc_None, "besselI0_diffraction", "besselI0_docstring", 0
    );

    double_slit_diffraction = PyUFunc_FromFuncAndData(
        double_slit_funcs, PyuFunc_None_3, five_real_in_one_real_out,
        3, 5, 1, PyUFunc_None, "double_slit_diffraction", 
        "double_slit_diffraction_docstring", 0
    );

    frequency_to_wavelength = PyUFunc_FromFuncAndData(
        frequency_to_wavelength_funcs, PyuFunc_None_3, one_real_in_one_real_out,
        3, 1, 1, PyUFunc_None, "frequency_to_wavelength",
        "frequency_to_wavelength_docstring", 0
    );

    fresnel_dpsi_dphi = PyUFunc_FromFuncAndData(
        dpsi_funcs, PyuFunc_None_3, seven_real_in_one_real_out, 3, 7, 1,
        PyUFunc_None, "fresnel_dpsi_dphi",  "fresnel_dpsi_dphi_docstring", 0
    );

    fresnel_d2psi_dphi2 = PyUFunc_FromFuncAndData(
        d2psi_funcs, PyuFunc_None_3, seven_real_in_one_real_out, 3, 7, 1,
        PyUFunc_None, "fresnel_d2psi_dphi2",  "fresnel_d2psi_dphi2_docstring", 0
    );

    fresnel_dpsi_dphi_ellipse = PyUFunc_FromFuncAndData(
        dpsi_ellipse_funcs, PyuFunc_None_3, nine_real_in_one_real_out, 3, 9, 1,
        PyUFunc_None, "fresnel_dpsi_dphi_ellipse", 
        "fresnel_dpsi_dphi_ellipse_docstring", 0
    );

    fresnel_cos = PyUFunc_FromFuncAndData(
        fresnel_cos_funcs, PyuFunc_None_3, one_real_in_one_real_out, 3, 1, 1,
        PyUFunc_None, "fresnel_cos", "fresnel_cos_docstring", 0
    );

    fresnel_psi = PyUFunc_FromFuncAndData(
        psi_funcs, PyuFunc_None_3, seven_real_in_one_real_out, 3, 7, 1,
        PyUFunc_None, "fresnel_psi",  "fresnel_psi_docstring", 0
    );

    fresnel_scale = PyUFunc_FromFuncAndData(
        fresnel_scale_funcs, PyuFunc_None_3, four_real_in_one_real_out, 3, 4, 1,
        PyUFunc_None, "fresnel_scale", "fresnel_scale_docstring", 0
    );

    fresnel_sin = PyUFunc_FromFuncAndData(
        fresnel_sin_funcs, PyuFunc_None_3, one_real_in_one_real_out, 3, 1, 1,
        PyUFunc_None, "fresnel_sin", "fresnel_sin_docstring", 0
    );

    lambertw = PyUFunc_FromFuncAndData(
        lambertw_funcs, PyuFunc_None_3, one_real_in_one_real_out, 3, 1, 1,
        PyUFunc_None, "lambertw", "lambertw_docstring", 0
    );

    resolution_inverse = PyUFunc_FromFuncAndData(
        res_inv_funcs, PyuFunc_None_3, one_real_in_one_real_out,
        3, 1, 1, PyUFunc_None, "resolution_inverse", 
        "resolution_inverse_docstring", 0
    );

    sinc = PyUFunc_FromFuncAndData(
        sinc_funcs, PyuFunc_None_3, one_real_in_one_real_out,
        3, 1, 1, PyUFunc_None, "sinc",  "sinc_docstring", 0
    );

    single_slit_diffraction = PyUFunc_FromFuncAndData(
        single_slit_funcs, PyuFunc_None_3, three_real_in_one_real_out,
        3, 3, 1, PyUFunc_None, "single_slit_diffraction", 
        "single_slit_diffraction_docstring", 0
    );

    wavelength_to_wavenumber = PyUFunc_FromFuncAndData(
        wavelength_to_wavenumber_funcs, PyuFunc_None_3,
        one_real_in_one_real_out, 3, 1, 1, PyUFunc_None,
        "wavelength_to_wavenumber", "wavelength_to_wavenumber_docstring", 0
    );

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "besselJ0", besselJ0);
    PyDict_SetItemString(d, "besselI0", besselI0);
    PyDict_SetItemString(d, "double_slit_diffraction", double_slit_diffraction);
    PyDict_SetItemString(d, "frequency_to_wavelength", frequency_to_wavelength);
    PyDict_SetItemString(d, "fresnel_cos", fresnel_cos);
    PyDict_SetItemString(d, "fresnel_psi", fresnel_psi);
    PyDict_SetItemString(d, "fresnel_dpsi_dphi", fresnel_dpsi_dphi);
    PyDict_SetItemString(d, "fresnel_d2psi_dphi2", fresnel_d2psi_dphi2);
    PyDict_SetItemString(d, "fresnel_dpsi_dphi_ellipse",
                         fresnel_dpsi_dphi_ellipse);
    PyDict_SetItemString(d, "fresnel_scale", fresnel_scale);
    PyDict_SetItemString(d, "fresnel_sin", fresnel_sin);
    PyDict_SetItemString(d, "lambertw", lambertw);
    PyDict_SetItemString(d, "resolution_inverse", resolution_inverse);
    PyDict_SetItemString(d, "sinc", sinc);
    PyDict_SetItemString(d, "single_slit_diffraction", single_slit_diffraction);
    PyDict_SetItemString(d, "wavelength_to_wavenumber",
                         wavelength_to_wavenumber);

    Py_DECREF(double_slit_diffraction);
    Py_DECREF(besselJ0);
    Py_DECREF(besselI0);
    Py_DECREF(frequency_to_wavelength);
    Py_DECREF(fresnel_cos);
    Py_DECREF(fresnel_psi);
    Py_DECREF(fresnel_dpsi_dphi);
    Py_DECREF(fresnel_d2psi_dphi2);
    Py_DECREF(fresnel_dpsi_dphi_ellipse);
    Py_DECREF(fresnel_scale);
    Py_DECREF(fresnel_sin);
    Py_DECREF(lambertw);
    Py_DECREF(resolution_inverse);
    Py_DECREF(sinc);
    Py_DECREF(single_slit_diffraction);
    Py_DECREF(wavelength_to_wavenumber);

    return m;
}
#endif