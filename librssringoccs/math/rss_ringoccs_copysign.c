#include <rss_ringoccs/include/rss_ringoccs_math.h>

double rssringoccs_Double_Copysign(double x, double y)
{
    double out;

    if (y < 0)
    {
        if (x < 0)
            out = x;
        else
            out = -x;
    }
    else if (0 < y)
    {
        if (x < 0)
            out = -x;
        else
            out = x;
    }
    else
        out = 0;

    return out;
}
