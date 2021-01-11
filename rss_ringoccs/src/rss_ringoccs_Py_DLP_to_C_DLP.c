
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
static double *__extract_data(rssringoccs_DLPObj *dlp, PyObject *py_dlp,
                              const char *var_name)
{
    PyObject *tmp;
    PyObject *arr;
    unsigned long len;

    if (dlp == NULL)
        return NULL;

    if (dlp->error_occurred)
        return NULL;

    if (!PyObject_HasAttrString(py_dlp, var_name))
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = malloc(sizeof(*dlp->error_message) * 256);
        sprintf(
            dlp->error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rInput DLP Instance is missing the following attribute:\n"
            "\r\t%s\n\n",
            var_name
        );
        return NULL;
    }
    else
        tmp = PyObject_GetAttrString(py_dlp, var_name);

    if (!PyArray_Check(tmp))
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = malloc(sizeof(*dlp->error_message) * 256);
        sprintf(
            dlp->error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r%s must be a numpy array.\n",
            var_name
        );
        return NULL;
    }
    else
        arr = PyArray_FromObject(tmp, NPY_DOUBLE, 1, 1);

    len = (unsigned long)PyArray_DIMS((PyArrayObject *)arr)[0];

    /*  If PyArray_FromObject failed arr should be NULL. If so, raise error.  */
    if (!arr)
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = malloc(sizeof(*dlp->error_message) * 256);
        sprintf(
            dlp->error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r%s must be a numpy array.\n",
            var_name
        );
        return NULL;
    }

    /*  Currently we only allow for one dimensional inputs.                   */
    else if (PyArray_NDIM((PyArrayObject *)arr) != 1)
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = malloc(sizeof(*dlp->error_message) * 256);
        sprintf(
            dlp->error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r%s must be a one-dimensional numpy array.\n",
            var_name
        );
        return NULL;
    }

    /*  arr should have the same number of elements as rho_km_vals.           */
    else if (len != dlp->arr_size)
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = malloc(sizeof(*dlp->error_message) * 256);
        sprintf(
            dlp->error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r%s and rho_km_vals have a different number of elements.\n",
            var_name
        );
        return NULL;
    }

    /*  If every passed, set ptr to point to the data inside the array arr.   */
    return (double *)PyArray_DATA((PyArrayObject *)arr);;
}

rssringoccs_DLPObj *rssringoccs_Py_DLP_to_C_DLP(PyObject *py_dlp)
{
    PyObject *tmp;
    PyObject *arr;
    rssringoccs_DLPObj *dlp;

    if (py_dlp == NULL)
        return NULL;

    dlp = malloc(sizeof(*dlp));
    if (dlp == NULL)
        return dlp;

    dlp->error_occurred = rssringoccs_False;
    dlp->error_message = NULL;

    if (py_dlp == NULL)
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rInput DLP Instance is NULL.\n"
        );
        return dlp;
    }

    /*  Next we're going to run error checks on the input numpy arrays which  *
     *  should be contained inside of the DLPInst object. We'll check that    *
     *  these attributes exist, that they are numpy arrays, are 1 dimensional,*
     *  and have the same number of elements as rho_km_vals. We'll also       *
     *  convert the arrays to double and retrieve a pointer to the data.      *
     *  First, we need to make sure rho_km_vals is a legal numpy array and    *
     *  extract the length of it. Check that rho_km_vals exists in DLPInst.   */
    if (!PyObject_HasAttrString(py_dlp, "rho_km_vals"))
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rInput DLP Instance is missing the following attribute:\n"
            "\r\trho_km_vals\n\n"
        );
        return dlp;
    }

    /*  If it exists, get a pointer to it.                                    */
    else
        tmp = PyObject_GetAttrString(py_dlp, "rho_km_vals");

    /*  Now make sure rho_km_vals is a numpy array.                           */
    if (!PyArray_Check(tmp))
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rrho_km_vals must be a numpy array.\n"
        );
        return dlp;
    }

    /*  If rho_km_vals is a numpy array, try to convert it to double.         */
    else
        arr = PyArray_FromObject(tmp, NPY_DOUBLE, 1, 1);

    /*  If PyArray_FromObject failed arr should be NULL. If so, raise error. */
    if (!arr)
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rCould not convert rho_km_vals to double array. Input is most\n"
            "\rlikely complex numbers or contains a string.\n\n"
        );
        return dlp;
    }

    /*  Currently we only allow for one dimensional inputs.                   */
    else if (PyArray_NDIM((PyArrayObject *)arr) != 1)
    {
        dlp->error_occurred = rssringoccs_True;
        dlp->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rrho_km_vals must be a one-dimensional numpy array.\n"
        );
        return dlp;
    }

    /*  If every passed, set tau.rho_km_vals to point to the data inside arr. */
    dlp->rho_km_vals = (double *)PyArray_DATA((PyArrayObject *)arr);
    dlp->arr_size = PyArray_DIMS((PyArrayObject *)arr)[0];

    dlp->p_norm_vals      = __extract_data(dlp, py_dlp, "p_norm_vals");
    dlp->phase_rad_vals   = __extract_data(dlp, py_dlp, "phase_rad_vals");
    dlp->phi_rad_vals     = __extract_data(dlp, py_dlp, "phi_rad_vals");
    dlp->phi_rl_rad_vals  = __extract_data(dlp, py_dlp, "phi_rl_rad_vals");
    dlp->B_rad_vals       = __extract_data(dlp, py_dlp, "B_rad_vals");
    dlp->D_km_vals        = __extract_data(dlp, py_dlp, "D_km_vals");
    dlp->f_sky_hz_vals    = __extract_data(dlp, py_dlp, "f_sky_hz_vals");
    dlp->rho_dot_kms_vals = __extract_data(dlp, py_dlp, "rho_dot_kms_vals");
    dlp->t_oet_spm_vals   = __extract_data(dlp, py_dlp, "t_oet_spm_vals");
    dlp->t_ret_spm_vals   = __extract_data(dlp, py_dlp, "t_ret_spm_vals");
    dlp->t_set_spm_vals   = __extract_data(dlp, py_dlp, "t_set_spm_vals");
    dlp->rx_km_vals       = __extract_data(dlp, py_dlp, "rx_km_vals");
    dlp->ry_km_vals       = __extract_data(dlp, py_dlp, "ry_km_vals");
    dlp->rz_km_vals       = __extract_data(dlp, py_dlp, "rz_km_vals");

    dlp->rho_corr_pole_km_vals
        = __extract_data(dlp, py_dlp, "rho_corr_pole_km_vals");

    dlp->rho_corr_timing_km_vals
        = __extract_data(dlp, py_dlp, "rho_corr_timing_km_vals");

    dlp->raw_tau_threshold_vals
        = __extract_data(dlp, py_dlp, "raw_tau_threshold_vals");

    return dlp;
}
