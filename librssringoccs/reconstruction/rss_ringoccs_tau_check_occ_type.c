#include <stdlib.h>
#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

void rssringoccs_Tau_Check_Occ_Type(rssringoccs_TAUObj *tau)
{
    double min_val;
    double max_val;
    unsigned long n;

    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    if (tau->rho_dot_kms_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Occ_Type\n\n"
            "\rInput rho_dot_kms_vals is NULL. Returning.\n"
        );
        return;
    }

    if (tau->arr_size <= 2)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Copy_DLP_Data_To_Tau\n\n"
            "\rInput arrays have less than 2 points. It is impossible to\n"
            "\rperform reconstrunction. Returning.\n"
        );
        return;
    }

    min_val = rssringoccs_Min_Double(tau->rho_dot_kms_vals, tau->arr_size);
    max_val = rssringoccs_Max_Double(tau->rho_dot_kms_vals, tau->arr_size);

    if ((min_val < 0.0) && (max_val > 0.0))
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Occ_Type\n\n"
            "\r\tdrho/dt has positive and negative values.\n"
            "\r\tYour input file is probably a chord occultation.\n"
            "\r\tDiffraction Correction can only be performed for\n"
            "\r\tone event at a time. That is, ingress or egress.\n\n"
            "\r\tTO CORRECT THIS:\n"
            "\r\t\tSplit the input into two parts: Ingress and Engress\n"
            "\r\t\tand perform diffraction correction twice.\n\n"
        );
        return;
    }
    else if ((min_val == 0.0) || (max_val == 0.0))
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Occ_Type\n\n"
            "\r\tdrho/dt has zero valued elements.\n"
            "\r\tYour input file is probably a chord occultation.\n"
            "\r\tDiffraction Correction can only be performed for\n"
            "\r\tone event at a time. That is, ingress or egress.\n\n"
            "\r\tTO CORRECT THIS:\n"
            "\r\t\tSplit the input into two parts: Ingress and Engress\n"
            "\r\t\tand perform diffraction correction twice.\n\n"
        );
        return;
    }
    else if ((tau->dx_km > 0.0) && (max_val < 0.0))
    {
        for(n=0; n<tau->arr_size; ++n)
            tau->rho_dot_kms_vals[n]
                = rssringoccs_Double_Abs(tau->rho_dot_kms_vals[n]);
    }
    else if ((tau->dx_km < 0.0) && (min_val > 0.0))
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_Occ_Type\n\n"
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r\trho_km_vals is decreasing yet rho_dot_kms_vals\n"
            "\r\tis positiive. Check DLP class for errors.\n\n"
        );
        return;
    }
    else if (tau->dx_km < 0.0)
    {
        rssringoccs_Reverse_Double_Array(tau->rho_km_vals,      tau->arr_size);
        rssringoccs_Reverse_Double_Array(tau->phi_rad_vals,     tau->arr_size);
        rssringoccs_Reverse_Double_Array(tau->B_rad_vals,       tau->arr_size);
        rssringoccs_Reverse_Double_Array(tau->D_km_vals,        tau->arr_size);
        rssringoccs_Reverse_Double_Array(tau->phase_rad_vals,   tau->arr_size);
        rssringoccs_Reverse_Double_Array(tau->p_norm_vals,      tau->arr_size);
        rssringoccs_Reverse_Double_Array(tau->f_sky_hz_vals,    tau->arr_size);
        rssringoccs_Reverse_Double_Array(tau->rho_dot_kms_vals, tau->arr_size);

        for(n=0; n<tau->arr_size; ++n)
            tau->rho_dot_kms_vals[n]
                = rssringoccs_Double_Abs(tau->rho_dot_kms_vals[n]);

        tau->dx_km = -tau->dx_km;
    }
}
