#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

rssringoccs_ThreeVector
rssringoccs_Inverse_Orthographic_Projection_Torus(rssringoccs_TwoVector P,
                                                  double r, double R)
{
    /*  Declare all necessary variables. C89 requires this at the top.        */
    double x, y, z, threshold;
    rssringoccs_ThreeVector out;

    /*  Extract the X and Y components from the point P.                      */
    x = rssringoccs_TwoVector_X(P);
    y = rssringoccs_TwoVector_Y(P);

    threshold = rssringoccs_Double_Sqrt(x*x + y*y) - R;
    threshold = r*r - threshold*threshold;

    if (threshold < 0.0)
    {
        out = rssringoccs_ThreeVector_Rect(rssringoccs_NaN,
                                           rssringoccs_NaN,
                                           rssringoccs_NaN);
    }
    else
    {
        z = rssringoccs_Double_Sqrt(threshold);
        out = rssringoccs_ThreeVector_Rect(x, y, z);
    }

    return out;
}
