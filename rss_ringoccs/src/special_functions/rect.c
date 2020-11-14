/*  The C Standard Library header for math functions and more found here.     */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

#define _define_rect_window(type, Type)                                        \
type rssringoccs_Rect_Window_##Type(type x, type W)                            \
{                                                                              \
    type abs_x, rect_x;                                                        \
                                                                               \
    abs_x = rssringoccs_Abs_##Type(x);                                         \
                                                                               \
    if (abs_x <= 0.5*W)                                                        \
        rect_x = 1.0;                                                          \
    else                                                                       \
        rect_x = 0.0;                                                          \
                                                                               \
    return rect_x;                                                             \
}

_define_rect_window(float, Float)
_define_rect_window(double, Double)
_define_rect_window(long double, LongDouble)

