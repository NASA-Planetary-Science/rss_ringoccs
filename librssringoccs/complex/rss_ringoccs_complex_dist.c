

#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

double
rssringoccs_CDouble_Dist(rssringoccs_ComplexDouble z0,
                         rssringoccs_ComplexDouble z1)
{
    double x0, y0, x1, y1, dx, dy, dist;

    x0 = rssringoccs_CDouble_Real_Part(z0);
    y0 = rssringoccs_CDouble_Imag_Part(z0);
    x1 = rssringoccs_CDouble_Real_Part(z1);
    y1 = rssringoccs_CDouble_Imag_Part(z1);

    dx = x1-x0;
    dy = y1-y0;

    dist = rssringoccs_Double_Sqrt(dx*dx + dy*dy);
    return dist;
}
