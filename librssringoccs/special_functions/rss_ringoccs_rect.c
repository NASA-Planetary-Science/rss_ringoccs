/*  The C Standard Library header for math functions and more found here.     */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

#define _define_rect_window(type, Type)                                        \
type rssringoccs_##Type##_Rect_Window(type x, type W)                          \
{                                                                              \
    type abs_x, rect_x;                                                        \
                                                                               \
    abs_x = rssringoccs_##Type##_Abs(x);                                       \
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
_define_rect_window(long double, LDouble)
