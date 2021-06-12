/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_RSR_READER_H
#define RSS_RINGOCCS_RSR_READER_H

#include <limits.h>

#if CHAR_BIT == 8
typedef char rssringoccs_int_8;
typedef signed char rssringoccs_sint_8;
typedef unsigned char rssringoccs_uint_8;
#else
#error "char is not 8 bits. Cannot read RSR files."
#endif

#if USHRT_MAX == 0xFFFF
typedef unsigned short int rssringoccs_uint_16;
#elif UINT_MAX == 0xFFFF
typedef unsigned int rssringoccs_uint_16;
#else
#error "neither short nor int is 16 bits. Cannot read RSR files."
#endif

#if UINT_MAX == 0xFFFFFFFF
typedef int rssringoccs_int_32;
typedef unsigned int rssringoccs_uint_32;
#elif ULONG_MAX == 0xFFFFFFFF
typedef int rssringoccs_int_32;
typedef unsigned long int rssringoccs_uint_32; 
#error "neither int nor long is 32 bits. Cannot read RSR files."
#endif

#if ULONG_MAX == 0xFFFFFFFFFFFFFFFF
typedef unsigned long int rssringoccs_uint_64;
#elif defined(ULONG_MAX) && ULONG_MAX == 0xFFFFFFFFFFFFFFFF
typedef unsigned long long int rssringoccs_uint_64;
#else 
#error "neither long nor long long is 64 bits. Cannot read RSR files."
#endif

typedef struct rssringoccs_SFDUObj {
    rssringoccs_int_8 authority[4];
    rssringoccs_int_8 version;
    rssringoccs_int_8 class;
    rssringoccs_int_8 reserved[2];
    rssringoccs_int_8 data_description[4];
    rssringoccs_uint_64 length;
} rssringoccs_SFDUObj;

typedef struct rssringoccs_HeaderAggregationObj {
    rssringoccs_uint_16 type;
    rssringoccs_uint_16 length;
} rssringoccs_HeaderAggregationObj;

typedef struct rssringoccs_RSRPrimaryHeaderObj {
    rssringoccs_uint_16 type;
    rssringoccs_uint_16 length;
    rssringoccs_uint_8 data_major;
    rssringoccs_uint_8 data_minor;
    rssringoccs_uint_8 mission_id;
    rssringoccs_uint_8 format_code;
} rssringoccs_RSRPrimaryHeaderObj;

typedef struct rssringoccs_RSRSecondaryHeaderObj {
    rssringoccs_int_8 reserved;
    rssringoccs_int_8 ul_band;
    rssringoccs_int_8 dl_band;
    rssringoccs_sint_8 fgain_px_no;
    rssringoccs_uint_8 originator;
    rssringoccs_uint_8 last_modifier;
    rssringoccs_uint_8 spc_id;
    rssringoccs_uint_8 dss_id;
    rssringoccs_uint_8 rsr_id;
    rssringoccs_uint_8 schan_id;
    rssringoccs_uint_8 spacecraft;
    rssringoccs_uint_8 trk_mode;
    rssringoccs_uint_8 ul_dss_id;
    rssringoccs_uint_8 fgain_if_bandwidth;
    rssringoccs_uint_8 frov_flag;
    rssringoccs_uint_8 attenuation;
    rssringoccs_uint_8 adc_rms;
    rssringoccs_uint_8 adc_peak;
    rssringoccs_uint_8 bits_per_sample;
    rssringoccs_uint_8 data_error;
    rssringoccs_uint_16 type;
    rssringoccs_uint_16 length;
    rssringoccs_uint_16 rsr_software_id;
    rssringoccs_uint_16 record_sequence_number;
    rssringoccs_uint_16 prdx_pass_number;
    rssringoccs_uint_16 year;
    rssringoccs_uint_16 doy;
    rssringoccs_uint_16 sample_rate;
    rssringoccs_uint_16 ddc_lo;
    rssringoccs_uint_16 rfif_lo;
    rssringoccs_uint_16 sfdu_year;
    rssringoccs_uint_16 sfdu_doy;
    rssringoccs_uint_32 seconds;
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
    rssringoccs_uint_16 type;
    rssringoccs_uint_16 length;
    rssringoccs_int_32 *data;
} rssringoccs_RSRDataObj;

typedef struct rssringoccs_RDRObj {
  rssringoccs_SFDUObj SFDUH;
  rssringoccs_HeaderAggregationObj HA;
  rssringoccs_RSRPrimaryHeaderObj PH;
  rssringoccs_RSRSecondaryHeaderObj SH;
  rssringoccs_RSRDataObj DATA;
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

extern rssringoccs_RDRObj *
rssringoccs_Get_RDR_Object(const char *rsr_file);

#endif
/*  End of include guard.                                                     */
