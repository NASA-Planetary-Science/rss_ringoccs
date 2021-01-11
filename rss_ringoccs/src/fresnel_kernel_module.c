


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
