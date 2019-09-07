/*  To avoid compiler warnings about deprecated numpy stuff.                  */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/*  complex data types, as well as _Complex_I, are defined here.              */
#include <complex.h>

/*  compute_norm_eq, max and min found here. math.h included here as well.    */
#include "__math_functions.h"

/* Include fresnel integrals header. This includes frensel_sin/cos.           */
#include "_fraunhofer_diffraction_wrappers.h"
#include "_fresnel_diffraction_wrappers.h"
#include "_fresnel_integrals_wrappers.h"
#include "_fresnel_kernel_wrappers.h"
#include "_lambertw_wrappers.h"
#include "_physics_functions_wrappers.h"
#include "_resolution_inverse_function_wrappers.h"
#include "_sinc_wrappers.h"
#include "_math_function_wrappers.h"

/*  Various header files required for the C-Python API to work.               */
#include <Python.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

/*---------------------------DEFINE PYTHON FUNCTIONS-------------------------*
 * This contains the Numpy-C and Python-C API parts that allow for the above *
 * functions to be called in Python. Numpy arrays, as well as floating point *
 * and integer valued arguments may then be passed into these functions for  *
 * improvement in performance, as opposed to the routines written purely in  *
 * Python. Successful compiling requires the Numpy and Python header files.  *
 *---------------------------------------------------------------------------*/
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
                         "\n\r\trss_ringoccs.diffrec.math_functions.max\n"
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
                         "\n\r\trss_ringoccs.diffrec.math_functions.max\n"
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
                         "\n\r\trss_ringoccs.diffrec.math_functions.max\n"
                         "\r\t\tInput must be one dimensional.");
            return NULL;
        }

        // Useful information about the data.
        typenum = PyArray_TYPE(arr);
        dim     = PyArray_DIMS(arr)[0];
        data    = PyArray_DATA(arr);

        if (dim == 0){
            PyErr_Format(PyExc_TypeError,
                         "\n\r\trss_ringoccs.diffrec.math_functions.max\n"
                         "\r\t\tInput is zero dimensional.");
            return NULL;
        }

        if (typenum == NPY_FLOAT){
            return PyFloat_FromDouble(
                Max_Float((float *)data, dim)
            );
        }
        else if (typenum == NPY_DOUBLE){
            return PyFloat_FromDouble(
                Max_Double((double *)data, dim)
            );
        }
        else if (typenum == NPY_LONGDOUBLE){
            return PyFloat_FromDouble(
                Max_Long_Double((long double *)data, dim)
            );
        }
        else if (typenum == NPY_SHORT){
            return PyLong_FromLong(
                Max_Short((short *)data, dim)
            );
        }
        else if (typenum == NPY_INT){
            return PyLong_FromLong(
                Max_Int((int *)data, dim)
            );
        }
        else if (typenum == NPY_LONG){
            return PyLong_FromLong(
                Max_Long((long *)data, dim)
            );
        }
        else if (typenum == NPY_LONGLONG){
            return PyLong_FromLong(
                Max_Long_Long((long long *)data, dim)
            );
        }
        else {
            PyErr_Format(PyExc_TypeError,
                         "\n\r\trss_ringoccs.diffrec.math_functions.max\n"
                         "\r\t\tInput should be a numpy array of numbers.");
            return NULL;
        }
    }
    else {
        PyErr_Format(PyExc_TypeError,
                     "\n\r\trss_ringoccs.diffrec.math_functions.max\n"
                     "\r\t\tInput should be a numpy array of numbers.");
        return NULL;
    }
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

static PyMethodDef _special_functions_methods[] =
{
    {"compute_norm_eq", compute_norm_eq,
     METH_VARARGS, "Compute the normalized equivalent width of an array."},
    {"max", max, METH_VARARGS, "Compute the maximum of a numpy array."},
    {"min", min, METH_VARARGS, "Compute the minimum of a numpy array."},
    {NULL, NULL, 0, NULL}
};
/*-------------------------DEFINE UNIVERSAL FUNCTIONS-------------------------*/
PyUFuncGenericFunction besselJ0_funcs[3] = {
    &float_besselJ0,
    &double_besselJ0,
    &long_double_besselJ0
};

PyUFuncGenericFunction besselI0_funcs[3] = {
    &float_besselI0,
    &double_besselI0,
    &long_double_besselI0
};

PyUFuncGenericFunction double_slit_funcs[3] = {
    &float_double_slit_diffraction,
    &double_double_slit_diffraction,
    &long_double_double_slit_diffraction
};

PyUFuncGenericFunction invsqwellsol_funcs[3] = {
    &complex_float_inv_square_well,
    &complex_double_inv_square_well,
    &complex_long_double_inv_square_well
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

PyUFuncGenericFunction left_straightedge_funcs[3] = {
    &float_left_straightedge,
    &double_left_straightedge,
    &long_double_left_straightedge
};

PyUFuncGenericFunction psi_funcs[3] = {
    &float_fresnel_psi,
    &double_fresnel_psi,
    &long_double_fresnel_psi
};

PyUFuncGenericFunction dpsi_funcs[3] = {
    float_fresnel_dpsi_dphi,
    &double_fresnel_dpsi_dphi,
    &long_double_fresnel_dpsi_dphi
};

PyUFuncGenericFunction res_inv_funcs[3] = {
    &float_resolution_inverse,
    &double_resolution_inverse,
    &long_double_resolution_inverse
};

PyUFuncGenericFunction right_straightedge_funcs[3] = {
    &float_right_straightedge,
    &double_right_straightedge,
    &long_double_right_straightedge
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

PyUFuncGenericFunction sqwellsol_funcs[3] = {
    &complex_float_square_well,
    &complex_double_square_well,
    &complex_long_double_square_well
};

PyUFuncGenericFunction sqwellphase_funcs[3] = {
    &float_square_well_phase,
    &double_square_well_phase,
    &long_double_square_well_phase
};

PyUFuncGenericFunction wavelength_to_wavenumber_funcs[3] = {
    &float_wavelength_to_wavenumber,
    &double_wavelength_to_wavenumber,
    &long_double_wavelength_to_wavenumber
};

/*  Input and return types for double input and out.                          */
static void *PyuFunc_None_3[3] = {NULL, NULL, NULL};

static char one_real_in_one_real_out[6] = {NPY_FLOAT, NPY_FLOAT,
                                           NPY_DOUBLE, NPY_DOUBLE,
                                           NPY_LONGDOUBLE, NPY_LONGDOUBLE};

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

static char four_real_in_one_complex_out[15] = {
    NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_CFLOAT,
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_CDOUBLE,
    NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE,
    NPY_CLONGDOUBLE
};

static char seven_real_in_one_real_out[24] = {
    NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
    NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
    NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE,
    NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE
};

static char three_real_in_one_complex_out[12] = {
    NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_CFLOAT,
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_CDOUBLE,
    NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE, NPY_CLONGDOUBLE
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
    PyObject *inverse_square_well_diffraction;
    PyObject *frequency_to_wavelength;
    PyObject *fresnel_cos;
    PyObject *fresnel_psi;
    PyObject *fresnel_scale;
    PyObject *fresnel_sin;
    PyObject *fresnel_dpsi_dphi;
    PyObject *lambertw;
    PyObject *left_straightedge;
    PyObject *resolution_inverse;
    PyObject *right_straightedge;
    PyObject *single_slit_diffraction;
    PyObject *sinc;
    PyObject *square_well_diffraction;
    PyObject *square_well_phase;
    PyObject *wavelength_to_wavenumber;
    PyObject *m, *d;

    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    import_array();
    import_umath();

    besselJ0 = PyUFunc_FromFuncAndData(
        besselJ0_funcs, PyuFunc_None_3, one_real_in_one_real_out,
        3, 1, 1, PyUFunc_None, "besselJ0_diffraction", 
        "besselJ0_docstring", 0
    );

    besselI0 = PyUFunc_FromFuncAndData(
        besselI0_funcs, PyuFunc_None_3, one_real_in_one_real_out,
        3, 1, 1, PyUFunc_None, "besselI0_diffraction", 
        "besselI0_docstring", 0
    );

    double_slit_diffraction = PyUFunc_FromFuncAndData(
        double_slit_funcs, PyuFunc_None_3, four_real_in_one_real_out,
        3, 4, 1, PyUFunc_None, "double_slit_diffraction", 
        "double_slit_diffraction_docstring", 0
    );

    inverse_square_well_diffraction = PyUFunc_FromFuncAndData(
        invsqwellsol_funcs, PyuFunc_None_3, four_real_in_one_complex_out,
        3, 4, 1, PyUFunc_None, "inverse_square_well_diffraction", 
        "inverse_square_well_diffraction_docstring", 0
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

    left_straightedge = PyUFunc_FromFuncAndData(
        left_straightedge_funcs, PyuFunc_None_3, three_real_in_one_complex_out,
        3, 3, 1, PyUFunc_None, "left_straightedge", 
        "left_straightedge_docstring", 0
    );

    resolution_inverse = PyUFunc_FromFuncAndData(
        res_inv_funcs, PyuFunc_None_3, one_real_in_one_real_out,
        3, 1, 1, PyUFunc_None, "resolution_inverse", 
        "resolution_inverse_docstring", 0
    );

    right_straightedge = PyUFunc_FromFuncAndData(
        right_straightedge_funcs, PyuFunc_None_3, three_real_in_one_complex_out,
        3, 3, 1, PyUFunc_None, "right_straightedge", 
        "right_straightedge_docstring", 0
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

    square_well_diffraction = PyUFunc_FromFuncAndData(
        sqwellsol_funcs, PyuFunc_None_3, four_real_in_one_complex_out,
        3, 4, 1, PyUFunc_None, "square_well_diffraction", 
        "square_well_diffraction_docstring", 0
    );

    square_well_phase = PyUFunc_FromFuncAndData(
        sqwellphase_funcs, PyuFunc_None_3, four_real_in_one_complex_out,
        3, 4, 1, PyUFunc_None, "square_well_phase", 
        "square_well_phase_docstring", 0
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
    PyDict_SetItemString(d, "inverse_square_well_diffraction",
                         inverse_square_well_diffraction);
    PyDict_SetItemString(d, "frequency_to_wavelength", frequency_to_wavelength);
    PyDict_SetItemString(d, "fresnel_cos", fresnel_cos);
    PyDict_SetItemString(d, "fresnel_psi", fresnel_psi);
    PyDict_SetItemString(d, "fresnel_dpsi_dphi", fresnel_dpsi_dphi);
    PyDict_SetItemString(d, "fresnel_scale", fresnel_scale);
    PyDict_SetItemString(d, "fresnel_sin", fresnel_sin);
    PyDict_SetItemString(d, "lambertw", lambertw);
    PyDict_SetItemString(d, "left_straightedge", left_straightedge);
    PyDict_SetItemString(d, "resolution_inverse", resolution_inverse);
    PyDict_SetItemString(d, "right_straightedge", right_straightedge);
    PyDict_SetItemString(d, "sinc", sinc);
    PyDict_SetItemString(d, "single_slit_diffraction", single_slit_diffraction);
    PyDict_SetItemString(d, "square_well_diffraction", square_well_diffraction);
    PyDict_SetItemString(d, "square_well_phase", square_well_phase);
    PyDict_SetItemString(d, "wavelength_to_wavenumber",
                         wavelength_to_wavenumber);

    Py_DECREF(double_slit_diffraction);
    Py_DECREF(besselJ0);
    Py_DECREF(besselI0);
    Py_DECREF(inverse_square_well_diffraction);
    Py_DECREF(frequency_to_wavelength);
    Py_DECREF(fresnel_cos);
    Py_DECREF(fresnel_psi);
    Py_DECREF(fresnel_dpsi_dphi);
    Py_DECREF(fresnel_scale);
    Py_DECREF(fresnel_sin);
    Py_DECREF(lambertw);
    Py_DECREF(left_straightedge);
    Py_DECREF(resolution_inverse);
    Py_DECREF(right_straightedge);
    Py_DECREF(sinc);
    Py_DECREF(single_slit_diffraction);
    Py_DECREF(square_well_diffraction);
    Py_DECREF(square_well_phase);
    Py_DECREF(wavelength_to_wavenumber);

    return m;
}
#else
PyMODINIT_FUNC init__funcs(void)
{
    PyObject *double_slit_diffraction;
    PyObject *inverse_square_well_diffraction;
    PyObject *frequency_to_wavelength;
    PyObject *fresnel_cos;
    PyObject *fresnel_psi;
    PyObject *fresnel_scale;
    PyObject *fresnel_sin;
    PyObject *fresnel_dpsi_dphi;
    PyObject *lambertw;
    PyObject *left_straightedge;
    PyObject *resolution_inverse;
    PyObject *right_straightedge;
    PyObject *single_slit_diffraction;
    PyObject *sinc;
    PyObject *square_well_diffraction;
    PyObject *square_well_phase;
    PyObject *wavelength_to_wavenumber;
    PyObject *m, *d;

    m = Py_InitModule("__funcs", _special_functions_methods);
    if (m == NULL) {
        return;
    }

    import_array();
    import_umath();

    double_slit_diffraction = PyUFunc_FromFuncAndData(
        double_slit_funcs, PyuFunc_None_3, four_real_in_one_real_out,
        3, 4, 1, PyUFunc_None, "double_slit_diffraction", 
        "double_slit_diffraction_docstring", 0
    );

    inverse_square_well_diffraction = PyUFunc_FromFuncAndData(
        invsqwellsol_funcs, PyuFunc_None_3, four_real_in_one_complex_out,
        3, 4, 1, PyUFunc_None, "inverse_square_well_diffraction", 
        "inverse_square_well_diffraction_docstring", 0
    );

    frequency_to_wavelength = PyUFunc_FromFuncAndData(
        frequency_to_wavelength_funcs, PyuFunc_None_3, one_real_in_one_real_out,
        3, 1, 1, PyUFunc_None, "frequency_to_wavelength",
        "frequency_to_wavelength_docstring", 0
    );

    fresnel_dpsi_dphi = PyUFunc_FromFuncAndData(
        dpsi_funcs, PyuFunc_data, octo_double_types, 1, 7, 1,
        PyUFunc_None, "fresnel_dpsi_dphi",  "fresnel_dpsi_dphi_docstring", 0
    );

    fresnel_cos = PyUFunc_FromFuncAndData(
        fresnel_cos_funcs, PyuFunc_None_3, one_real_in_one_real_out, 3, 1, 1,
        PyUFunc_None, "fresnel_cos", "fresnel_cos_docstring", 0
    );

    fresnel_psi = PyUFunc_FromFuncAndData(
        psi_funcs, PyuFunc_data, octo_double_types, 1, 7, 1,
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

    left_straightedge = PyUFunc_FromFuncAndData(
        left_straightedge_funcs, PyuFunc_None_3, three_real_in_one_complex_out,
        3, 3, 1, PyUFunc_None, "left_straightedge", 
        "left_straightedge_docstring", 0
    );

    resolution_inverse = PyUFunc_FromFuncAndData(
        res_inv_funcs, PyuFunc_None_3, one_real_in_one_real_out,
        3, 1, 1, PyUFunc_None, "resolution_inverse", 
        "resolution_inverse_docstring", 0
    );

    right_straightedge = PyUFunc_FromFuncAndData(
        right_straightedge_funcs, PyuFunc_None_3, three_real_in_one_complex_out,
        3, 3, 1, PyUFunc_None, "right_straightedge", 
        "right_straightedge_docstring", 0
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

    square_well_diffraction = PyUFunc_FromFuncAndData(
        sqwellsol_funcs, PyuFunc_None_3, four_real_in_one_complex_out,
        3, 4, 1, PyUFunc_None, "square_well_diffraction", 
        "square_well_diffraction_docstring", 0
    );

    square_well_phase = PyUFunc_FromFuncAndData(
        sqwellphase_funcs, PyuFunc_None_3, four_real_in_one_complex_out,
        3, 4, 1, PyUFunc_None, "square_well_phase", 
        "square_well_phase_docstring", 0
    );

    wavelength_to_wavenumber = PyUFunc_FromFuncAndData(
        wavelength_to_wavenumber_funcs, PyuFunc_None_3,
        one_real_in_one_real_out, 3, 1, 1, PyUFunc_None,
        "wavelength_to_wavenumber", "wavelength_to_wavenumber_docstring", 0
    );

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "double_slit_diffraction", double_slit_diffraction);
    PyDict_SetItemString(d, "inverse_square_well_diffraction",
                         inverse_square_well_diffraction);
    PyDict_SetItemString(d, "frequency_to_wavelength", frequency_to_wavelength);
    PyDict_SetItemString(d, "fresnel_cos", fresnel_cos);
    PyDict_SetItemString(d, "fresnel_psi", fresnel_psi);
    PyDict_SetItemString(d, "fresnel_dpsi_dphi", fresnel_dpsi_dphi);
    PyDict_SetItemString(d, "fresnel_scale", fresnel_scale);
    PyDict_SetItemString(d, "fresnel_sin", fresnel_sin);
    PyDict_SetItemString(d, "lambertw", lambertw);
    PyDict_SetItemString(d, "left_straightedge", left_straightedge);
    PyDict_SetItemString(d, "resolution_inverse", resolution_inverse);
    PyDict_SetItemString(d, "right_straightedge", right_straightedge);
    PyDict_SetItemString(d, "sinc", sinc);
    PyDict_SetItemString(d, "single_slit_diffraction", single_slit_diffraction);
    PyDict_SetItemString(d, "square_well_diffraction", square_well_diffraction);
    PyDict_SetItemString(d, "square_well_phase", square_well_phase);
    PyDict_SetItemString(d, "wavelength_to_wavenumber",
                         wavelength_to_wavenumber);

    Py_DECREF(double_slit_diffraction);
    Py_DECREF(inverse_square_well_diffraction);
    Py_DECREF(frequency_to_wavelength);
    Py_DECREF(fresnel_cos);
    Py_DECREF(fresnel_psi);
    Py_DECREF(fresnel_dpsi_dphi);
    Py_DECREF(fresnel_scale);
    Py_DECREF(fresnel_sin);
    Py_DECREF(lambertw);
    Py_DECREF(left_straightedge);
    Py_DECREF(resolution_inverse);
    Py_DECREF(right_straightedge);
    Py_DECREF(sinc);
    Py_DECREF(single_slit_diffraction);
    Py_DECREF(square_well_diffraction);
    Py_DECREF(square_well_phase);
    Py_DECREF(wavelength_to_wavenumber);

    return m;
}
#endif