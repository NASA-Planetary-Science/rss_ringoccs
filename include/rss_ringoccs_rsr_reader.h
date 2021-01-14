/*  Include guard to avoid importing this file twice.                         */
#ifndef __RSS_RINGOCCS_RSR_READER_H__
#define __RSS_RINGOCCS_RSR_READER_H__

#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

typedef struct rssringoccs_SFDUObj {
    char authority[4];
    char version;
    char class;
    char reserved[2];
    char data_description[4];
    unsigned long length;
} rssringoccs_SFDUObj;

typedef struct rssringoccs_HeaderAggregationObj {
    unsigned short type;
    unsigned short length;
} rssringoccs_HeaderAggregationObj;

typedef struct rssringoccs_RSRPrimaryHeaderObj {
    unsigned short type;
    unsigned short length;
    unsigned char data_major;
    unsigned char data_minor;
    unsigned char mission_id;
    unsigned char format_code;
} rssringoccs_RSRPrimaryHeaderObj;

typedef struct rssringoccs_RSRSecondaryHeaderObj {
    unsigned int seconds;
    unsigned char originator;
    unsigned char last_modifier;
    unsigned char spc_id;
    unsigned char dss_id;
    unsigned char rsr_id;
    unsigned char schan_id;
    unsigned char spacecraft;
    unsigned char trk_mode;
    unsigned char ul_dss_id;
    unsigned char fgain_if_bandwidth;
    unsigned char frov_flag;
    unsigned char attenuation;
    unsigned char adc_rms;
    unsigned char adc_peak;
    unsigned char bits_per_sample;
    unsigned char data_error;
    unsigned short type;
    unsigned short length;
    unsigned short rsr_software_id;
    unsigned short record_sequence_number;
    unsigned short prdx_pass_number;
    unsigned short year;
    unsigned short doy;
    unsigned short sample_rate;
    unsigned short ddc_lo;
    unsigned short rfif_lo;
    unsigned short sfdu_year;
    unsigned short sfdu_doy;
    signed char fgain_px_no;
    char reserved;
    char ul_band;
    char dl_band;
    double sfdu_seconds;
    double predicts_time_shift;
    double predicts_freq_override;
    double predicts_freq_rate;
    double predicts_freq_offset;
    double sub_channel_freq;
    double rf_freq_point_1;
    double rf_freq_point_2;
    double rf_freq_point_3;
    double schan_freq_point_1;
    double schan_freq_point_2;
    double schan_freq_point_3;
    double schan_freq_poly_coef_1;
    double schan_freq_poly_coef_2;
    double schan_freq_poly_coef_3;
    double schan_accum_phase;
    double schan_phase_poly_coef_1;
    double schan_phase_poly_coef_2;
    double schan_phase_poly_coef_3;
    double schan_phase_poly_coef_4;
    double reserved2[2];
} rssringoccs_RSRSecondaryHeaderObj;

typedef struct rssringoccs_RSRDataObj {
    unsigned short type;
    unsigned short length;
    int *data;
} rssringoccs_RSRDataObj;

typedef struct rssringoccs_RDRObj {
  struct rssringoccs_SFDUObj SFDUH;
  struct rssringoccs_HeaderAggregationObj HA;
  struct rssringoccs_RSRPrimaryHeaderObj PH;
  struct rssringoccs_RSRSecondaryHeaderObj SH;
  struct rssringoccs_RSRDataObj DATA;
} rssringoccs_RDRObj;

typedef struct rssringoccs_RSRHeaderObj {

    /*  Second past midnight values.                                          */
    double *spm_vals;

    /*  Day of year for the event.                                            */
    unsigned int doy;

    /*  Year of the event.                                                    */
    unsigned int year;

    /*  The DSN station.                                                      */
    char *dsn;

    /*  Band of the RSR file (S, X, or Ka).                                   */
    char *band;

    /*  The sample rate from the RSR file.                                    */
    double sample_rate_hz;

} rssringoccs_RSRHeaderObj;
/*  End of rssringoccs_RSRHeaderObj definition.                               */

/*  Function for extracting the data from an rsr_file, given its path.        */
extern rssringoccs_RSRHeaderObj *
rssringoccs_Get_RSR_Header(const char *rsr_file);

#endif
/*  End of include guard.                                                     */
