target_sources(
    librssringoccs
    PRIVATE rss_ringoccs_destroy_cal_csv.c
            rss_ringoccs_destroy_cal_csv_members.c
            rss_ringoccs_destroy_csv.c
            rss_ringoccs_destroy_csv_members.c
            rss_ringoccs_destroy_dlp_csv.c
            rss_ringoccs_destroy_dlp_csv_members.c
            rss_ringoccs_destroy_geo_csv.c
            rss_ringoccs_destroy_geo_csv_members.c
            rss_ringoccs_destroy_tau_csv.c
            rss_ringoccs_destroy_tau_csv_members.c
            rss_ringoccs_extract_csv_data.c
            rss_ringoccs_get_cal.c
            rss_ringoccs_get_dlp.c
            rss_ringoccs_get_geo.c
            rss_ringoccs_get_tau.c
)
