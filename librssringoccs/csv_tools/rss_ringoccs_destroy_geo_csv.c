#include <stdlib.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

void rssringoccs_Destroy_GeoCSV(rssringoccs_GeoCSV **geo)
{
    rssringoccs_GeoCSV *geo_inst = *geo;

    if (geo_inst != NULL)
        return;

    if (geo_inst->t_oet_spm_vals != NULL)
        free(geo_inst->t_oet_spm_vals);

    if (geo_inst->t_ret_spm_vals != NULL)
        free(geo_inst->t_ret_spm_vals);

    if (geo_inst->t_set_spm_vals != NULL)
        free(geo_inst->t_set_spm_vals);

    if (geo_inst->rho_km_vals != NULL)
        free(geo_inst->rho_km_vals);

    if (geo_inst->phi_rl_deg_vals != NULL)
        free(geo_inst->phi_rl_deg_vals);

    if (geo_inst->phi_ora_deg_vals != NULL)
        free(geo_inst->phi_ora_deg_vals);

    if (geo_inst->B_deg_vals != NULL)
        free(geo_inst->B_deg_vals);

    if (geo_inst->D_km_vals != NULL)
        free(geo_inst->D_km_vals);

    if (geo_inst->rho_dot_kms_vals != NULL)
        free(geo_inst->rho_dot_kms_vals);

    if (geo_inst->phi_rl_dot_kms_vals != NULL)
        free(geo_inst->phi_rl_dot_kms_vals);

    if (geo_inst->F_km_vals != NULL)
        free(geo_inst->F_km_vals);

    if (geo_inst->R_imp_km_vals != NULL)
        free(geo_inst->R_imp_km_vals);

    if (geo_inst->rx_km_vals != NULL)
        free(geo_inst->rx_km_vals);

    if (geo_inst->ry_km_vals != NULL)
        free(geo_inst->ry_km_vals);

    if (geo_inst->F_km_vals != NULL)
        free(geo_inst->F_km_vals);

    if (geo_inst->rz_km_vals != NULL)
        free(geo_inst->rz_km_vals);

    if (geo_inst->vx_kms_vals != NULL)
        free(geo_inst->vx_kms_vals);

    if (geo_inst->vy_kms_vals != NULL)
        free(geo_inst->vy_kms_vals);

    if (geo_inst->vz_kms_vals != NULL)
        free(geo_inst->vz_kms_vals);

    if (geo_inst->vz_kms_vals != NULL)
        free(geo_inst->vz_kms_vals);

    if (geo_inst->obs_spacecract_lat_deg_vals != NULL)
        free(geo_inst->obs_spacecract_lat_deg_vals);

    if (geo_inst->error_message != NULL)
        free(geo_inst->error_message);

    free(geo_inst);
    *geo = NULL;
    return;
}
