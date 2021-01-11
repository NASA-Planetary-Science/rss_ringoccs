/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

#define _define_normeq(intype, InType, outtype)                                \
outtype rssringoccs_Normeq_##InType(intype *w_func, long N_elements)           \
{                                                                              \
    outtype out;                                                               \
    outtype tot_sq = 0.0;                                                      \
    outtype sq_tot = 0.0;                                                      \
    long n;                                                                    \
                                                                               \
    /* Compute both the sum and the sum of the squares.                      */\
    for (n = 0; n < N_elements; n++)                                           \
    {                                                                          \
        tot_sq += w_func[n];                                                   \
        sq_tot += w_func[n]*w_func[n];                                         \
    }                                                                          \
                                                                               \
    /*  Square the sum.                                                      */\
    tot_sq *= tot_sq;                                                          \
                                                                               \
    /*  Compute the normalized equivalent width and return.                  */\
    out =  N_elements * sq_tot / tot_sq;                                       \
    return out;                                                                \
}

_define_normeq(float, Float, float)
_define_normeq(double, Double, double)
_define_normeq(long double, LDouble, long double)
_define_normeq(int, Int, double)
_define_normeq(short, Short, double)
_define_normeq(long, Long, double)
