

void rssringoccs_Get_Py_Vars_From_Self(rssringoccs_TAUObj *tau,
                                       PyDiffrecObj *self)
{
    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    if (self == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_C_Tau_to_Py_Tau\n\n"
            "\rInput self is NULL. Aborting.n"
        );
        return;
    }

    tau->sigma    = self->sigma;
    tau->bfac     = self->bfac;
    tau->ecc      = self->ecc;
    tau->peri     = self->peri;
    tau->interp   = self->interp;
    tau->use_fwd  = self->use_fwd;
    tau->use_norm = self->use_norm;
    tau->verbose  = self->verbose;
}