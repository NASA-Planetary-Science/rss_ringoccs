/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_CALIBRATION_H
#define RSS_RINGOCCS_CALIBRATION_H

#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_complex.h>

/*  size_t typedef is given here.                                             */
#include <stdlib.h>

/*  Structure that contains all of the necessary data.                        */
typedef struct rssringoccs_DLPObj_Def {
    double *rho_km_vals;
    double *phi_deg_vals;
    double *B_deg_vals;
    double *D_km_vals;
    double *f_sky_hz_vals;
    double *rho_dot_kms_vals;
    double *t_oet_spm_vals;
    double *t_ret_spm_vals;
    double *t_set_spm_vals;
    double *rho_corr_pole_km_vals;
    double *rho_corr_timing_km_vals;
    double *phi_rl_deg_vals;
    double *p_norm_vals;
    double *phase_deg_vals;
    double *raw_tau_threshold_vals;
    double *rx_km_vals;
    double *ry_km_vals;
    double *rz_km_vals;
    size_t arr_size;
    tmpl_Bool error_occurred;
    char *error_message;
} rssringoccs_DLPObj;

typedef struct rssringoccs_CalcFreqOffsetObj_Def {

    /*  Observed event time at full sampling.                                 */
    double *t_oet_spm_vals;

    /*  Uncorrected real and imaginary components of signal.                  */
    tmpl_ComplexDouble *IQ_m;

    /*  Raw time sampling from spm_vals.                                      */
    double dt;

    /*  Half the width of the FFT window.                                     */
    double dt_freq;

    /*  Minimum time for sampling.                                            */
    double t_oet_spm_min;

    /*  Maximum time for sampling.                                            */
    double t_oet_spm_max;

    /*  Observed event time for frequency.                                    */
    double *f_oet_spm_vals;

    /*  Frequency offset, or frequency at max power.                          */
    double *f_offset_hz_vals;
} rssringoccs_CalcFreqOffsetObj;
/*  End of rssringoccs_CalcFreqOffsetObj definition.                          */

#endif
