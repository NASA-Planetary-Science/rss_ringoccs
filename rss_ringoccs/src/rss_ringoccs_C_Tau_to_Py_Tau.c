
/*  This function frees the memory allocated to a pointer by malloc when the  *
 *  corresponding variable is destroyed at the Python level. Without this you *
 *  will have serious memory leaks, so do not remove!                         */
static void capsule_cleanup(PyObject *capsule)
{
    void *memory = PyCapsule_GetPointer(capsule, NULL);
    free(memory);
}

static void __set_var(PyObject **py_ptr, double **ptr, unsigned long len)
{
    PyObject *arr;
    PyObject *capsule;
    PyObject *tmp;
    long pylength = (long)len;

    arr     = PyArray_SimpleNewFromData(1, &pylength, NPY_DOUBLE, *ptr);
    capsule = PyCapsule_New((void *) (*ptr), NULL, capsule_cleanup);

    PyArray_SetBaseObject((PyArrayObject *)arr, capsule);

    tmp = *py_ptr;
    Py_INCREF(arr);
    *py_ptr = arr;
    Py_XDECREF(tmp);
}

static void __set_cvar(PyObject **py_ptr, rssringoccs_ComplexDouble **ptr,
                       unsigned long len)
{
    PyObject *arr;
    PyObject *capsule;
    PyObject *tmp;
    long pylength = (long)len;

    arr     = PyArray_SimpleNewFromData(1, &pylength, NPY_CDOUBLE, *ptr);
    capsule = PyCapsule_New((void *)(*ptr), NULL, capsule_cleanup);

    PyArray_SetBaseObject((PyArrayObject *)arr, capsule);

    tmp = *py_ptr;
    Py_INCREF(arr);
    *py_ptr = arr;
    Py_XDECREF(tmp);
}

static void rssringoccs_C_Tau_to_Py_Tau(PyDiffrecObj *py_tau,
                                        rssringoccs_TAUObj *tau)
{
    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    if (py_tau == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_C_Tau_to_Py_Tau\n\n"
            "\rInput py_tau is NULL. Aborting.n"
        );
        return;
    }

    __set_cvar(&py_tau->T_hat_vals,      &tau->T_in,             tau->arr_size);
    __set_cvar(&py_tau->T_vals,          &tau->T_out,            tau->arr_size);

    __set_var(&py_tau->rho_km_vals,      &tau->rho_km_vals,      tau->arr_size);
    __set_var(&py_tau->B_rad_vals,       &tau->B_rad_vals,       tau->arr_size);
    __set_var(&py_tau->D_km_vals,        &tau->D_km_vals,        tau->arr_size);
    __set_var(&py_tau->F_km_vals,        &tau->F_km_vals,        tau->arr_size);
    __set_var(&py_tau->f_sky_hz_vals,    &tau->f_sky_hz_vals,    tau->arr_size);
    __set_var(&py_tau->p_norm_vals,      &tau->p_norm_vals,      tau->arr_size);
    __set_var(&py_tau->phase_rad_vals,   &tau->phase_rad_vals,   tau->arr_size);
    __set_var(&py_tau->phase_vals,       &tau->phase_vals,       tau->arr_size);
    __set_var(&py_tau->phi_rad_vals,     &tau->phi_rad_vals,     tau->arr_size);
    __set_var(&py_tau->phi_rl_rad_vals,  &tau->phi_rl_rad_vals,  tau->arr_size);
    __set_var(&py_tau->power_vals,       &tau->power_vals,       tau->arr_size);
    __set_var(&py_tau->rho_dot_kms_vals, &tau->rho_dot_kms_vals, tau->arr_size);
    __set_var(&py_tau->t_oet_spm_vals,   &tau->t_oet_spm_vals,   tau->arr_size);
    __set_var(&py_tau->t_ret_spm_vals,   &tau->t_ret_spm_vals,   tau->arr_size);
    __set_var(&py_tau->t_set_spm_vals,   &tau->t_set_spm_vals,   tau->arr_size);
    __set_var(&py_tau->tau_vals,         &tau->tau_vals,         tau->arr_size);
    __set_var(&py_tau->w_km_vals,        &tau->w_km_vals,        tau->arr_size);
    __set_var(&py_tau->rx_km_vals,       &tau->rx_km_vals,       tau->arr_size);
    __set_var(&py_tau->ry_km_vals,       &tau->ry_km_vals,       tau->arr_size);
    __set_var(&py_tau->rz_km_vals,       &tau->rz_km_vals,       tau->arr_size);

    __set_var(&py_tau->raw_tau_threshold_vals,
              &tau->raw_tau_threshold_vals, tau->arr_size);

    __set_var(&py_tau->rho_corr_pole_km_vals,
              &tau->rho_corr_pole_km_vals, tau->arr_size);

    __set_var(&py_tau->rho_corr_timing_km_vals,
              &tau->rho_corr_timing_km_vals, tau->arr_size);

    __set_var(&py_tau->tau_threshold_vals,
              &tau->tau_threshold_vals, tau->arr_size);
}

/*  To edit:
    PyObject          *p_norm_fwd_vals;
    PyObject          *phase_fwd_vals;
    PyObject          *T_hat_fwd_vals;
    PyObject          *rev_info;
    PyObject          *dathist;
    PyObject          *history;
*/
