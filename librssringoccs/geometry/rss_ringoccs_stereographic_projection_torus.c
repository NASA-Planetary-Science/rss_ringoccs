#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

rssringoccs_TwoVector
rssringoccs_Stereographic_Projection_Torus(double u, double v,
                                           double r, double R)
{
    double cos_u, sin_u, cos_v, sin_v, t, x, y, Px, Py, Qx, Qy;
    rssringoccs_TwoVector out;

    sin_u = rssringoccs_Double_Sin(u);

    if (sin_u == 1.0)
    {
        x = rssringoccs_Infinity;
        y = rssringoccs_Infinity;
    }
    else
    {
        t = 1.0 / (1.0-sin_u);
        cos_u = rssringoccs_Double_Cos(u);
        cos_v = rssringoccs_Double_Cos(v);
        sin_v = rssringoccs_Double_Sin(v);

        Qx = R*cos_v;
        Qy = R*sin_v;

        Px = cos_v*(R + r*cos_u);
        Py = sin_v*(R + r*cos_u);

        x = t*Px + (1.0 - t)*Qx;
        y = t*Py + (1.0 - t)*Qy;
    }

    out = rssringoccs_TwoVector_Rect(x, y);
    return out;
}
