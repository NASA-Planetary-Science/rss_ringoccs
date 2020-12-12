#include <stdio.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_ppm_plot.h>

double coeffsa[32] = {
    0.0,  341.0,    0.0, 0.0, 0.0, 0.0,  67518.0,   0.0, 0.0, 0.0,
    0.0, -398505.0, 0.0, 0.0, 0.0, 0.0, -1060200.0, 0.0, 0.0, 0.0,
    0.0,  326895.0, 0.0, 0.0, 0.0, 0.0,  10602.0,   0.0, 0.0, 0.0,
    0.0, -19.0
};

double coeffsb[31] = {
    -19.0,     0.0, 0.0, 0.0, 0.0, -10602.0,   0.0, 0.0, 0.0, 0.0,
     326895.0, 0.0, 0.0, 0.0, 0.0,  1060200.0, 0.0, 0.0, 0.0, 0.0,
    -398505.0, 0.0, 0.0, 0.0, 0.0, -67518.0,   0.0, 0.0, 0.0, 0.0,
     341.0
};

static rssringoccs_ComplexDouble func(rssringoccs_ComplexDouble z)
{
    rssringoccs_ComplexDouble f, g, out;

    f = rssringoccs_CDouble_Poly_Real_Coeffs(coeffsa, 31, z);
    g = rssringoccs_CDouble_Poly_Real_Coeffs(coeffsb, 30, z);
    out = rssringoccs_CDouble_Divide(f, g);
    return out;
}

int main(void)
{
    FILE *fp;
    unsigned int x, y, n;
    double z_x, z_y, norm;
    double rcp_factor;
    rssringoccs_ComplexDouble z_0, z_1, z_2, dist;

    unsigned int size = 4*1024;
    unsigned int maxIterations = 12;

    const double x_min = -100.0;
    const double x_max =  100.0;
    const double y_min = -100.0;
    const double y_max =  100.0;
    const double eps = 1.0e-4;

    fp = fopen("julia_set.ppm", "w");
    fprintf(fp, "P6\n%d %d\n255\n", size, size);

    rcp_factor = 1.0/(size-1.0);

    for (y=0; y<size; ++y)
    {
        z_y = (size - y - 1.0) * (y_max - y_min) * rcp_factor + y_min;
        for (x=0; x<size; ++x)
        {
            z_x = x * (x_max - x_min) * rcp_factor + x_min;
            z_0 = rssringoccs_CDouble_Rect(z_x, z_y);

            for(n = 0; n < maxIterations; n++)
            {
                z_1 = func(z_0);
                z_2 = func(z_1);
                dist = rssringoccs_CDouble_Subtract(z_0, z_2);
                norm = rssringoccs_CDouble_Abs(dist);

                if(norm < eps)
                    break;
                else
                    z_0 = z_1;
            }

            if(n == maxIterations)
                rssringoccs_Color(0, 0, 0, fp);
            else
                rssringoccs_RGB_Linear_Gradient(n, 0, maxIterations, fp);
        }
        printf("%d %d\n", y, n);
    }
    return 0;
}
