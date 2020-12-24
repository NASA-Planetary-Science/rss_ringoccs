#include <rss_ringoccs/include/rss_ringoccs_geometry.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_ppm_plot.h>
#include <stdlib.h>
#include <stdio.h>

/*  Initials

Binary                      Decimal

111100111110000000          249728
100100001000000000          147968
101000001001000001          164417
110000001001100011          197219
101000001001010101          164437
100100001001001001          148041
100010001001000001          139841
001000001001000001           33345
001111111000000000           65024

*/

static double random_unit_interval(void)
{
    return (double)rand()/RAND_MAX;
}

int G[9] = {
    65024, 33345, 139841, 148041, 164437, 197219, 164417, 147968, 249728
};

static int Tracer(rssringoccs_ThreeVector o, rssringoccs_ThreeVector d,
                  double *t, rssringoccs_ThreeVector *n)
{
    int k, j, m;
    double s, b, c, q, p, o_z, d_z;
    rssringoccs_ThreeVector P, temp;

    o_z = rssringoccs_ThreeVector_Z(o);
    d_z = rssringoccs_ThreeVector_Z(d);

    *t = 1.0e9;
    m = 0;
    p = - o_z/d_z;

    if (p > 0.01)
    {
        *t = p;
        *n = rssringoccs_ThreeVector_Rect(0.0, 0.0, 1.0);
        m = 1;
    }

    for (k=0; k<18; ++k)
    {
        for (j=0; j<9; ++j)
        {
            if (G[j] & (1 << k))
            {
                temp = rssringoccs_ThreeVector_Rect(-k, 0.0, -j-4);
                P = rssringoccs_ThreeVector_Add(o, temp);
                b = rssringoccs_ThreeVector_Dot_Product(P, d);
                c = rssringoccs_ThreeVector_Dot_Product(P, P) - 1.0;
                q = b*b - c;

                if (q>0.0)
                {
                    s = -b - rssringoccs_Double_Sqrt(q);

                    if ((s<*t) && (s > 0.01))
                    {
                        *t = s;
                        temp = rssringoccs_ThreeVector_Scale(*t, d);
                        temp = rssringoccs_ThreeVector_Add(P, temp);
                        *n = rssringoccs_ThreeVector_Normalize(temp);
                        m = 2;
                    }
                }
            }
        }
    }
    return m;
}

static rssringoccs_ThreeVector Shader(rssringoccs_ThreeVector o,
                                      rssringoccs_ThreeVector d)
{
    double t, b, p, factor, d_z, h_x, h_y;
    rssringoccs_ThreeVector n, out, h, temp, l, r, shade;
    int m;

    /*  Initialize variables.                                                 */
    t = 0.0;
    n = rssringoccs_ThreeVector_Rect(0.0, 0.0, 0.0);

    m = Tracer(o, d, &t, &n);
    d_z = rssringoccs_ThreeVector_Z(d);

    if (m == 0)
    {
        factor = d_z*d_z;
        factor = factor*factor;
        out = rssringoccs_ThreeVector_Rect(1.0-factor, 0.7, 1.0);
        return out;
    }

    temp = rssringoccs_ThreeVector_Scale(t, d);
    h = rssringoccs_ThreeVector_Add(o, temp);

    temp = rssringoccs_ThreeVector_Rect(9.0 + 5*random_unit_interval(),
                                        9.0 + 5*random_unit_interval(), 46.0);

    temp = rssringoccs_ThreeVector_Add(temp,
                                       rssringoccs_ThreeVector_Scale(-1.0, h));

    l = rssringoccs_ThreeVector_Normalize(temp);

    factor = -2.0*rssringoccs_ThreeVector_Dot_Product(n, d);
    temp = rssringoccs_ThreeVector_Scale(factor, n);
    r = rssringoccs_ThreeVector_Add(d, temp);

    b = rssringoccs_ThreeVector_Dot_Product(l, n);

    if ((b<0.0) || Tracer(h, l, &t, &n))
        b = 0.0;

    if (b>0)
        p = pow(rssringoccs_ThreeVector_Dot_Product(l, r), 99.0);
    else
        p = 0.0;

    if (m == 1)
    {
        h = rssringoccs_ThreeVector_Scale(0.2, h);
        h_x = ceil(rssringoccs_ThreeVector_X(h));
        h_y = ceil(rssringoccs_ThreeVector_Y(h));

        if (((int)h_x + (int)h_y) & 1)
        {
            temp = rssringoccs_ThreeVector_Rect(3.0, 1.0, 1.0);
            factor = 0.2*b + 0.1;
            out = rssringoccs_ThreeVector_Scale(factor, temp);
            return out;
        }
        else
        {
            temp = rssringoccs_ThreeVector_Rect(3.0, 3.0, 3.0);
            factor = 0.2*b + 0.1;
            out = rssringoccs_ThreeVector_Scale(factor, temp);
            return out;
        }
    }

    temp = rssringoccs_ThreeVector_Rect(p, p, p);
    shade = rssringoccs_ThreeVector_Scale(0.5, Shader(h, r));

    out = rssringoccs_ThreeVector_Add(shade, temp);
    return out;
}

int main(void)
{
    rssringoccs_ThreeVector g, a, b, c, p, t, temp1, temp2, alpha, beta;
    unsigned int x, y, r;
    double factor, px, py, pz;
    FILE *fp;
    unsigned int size = 1024;

    g = rssringoccs_ThreeVector_Rect(-6.0, -16.0, 0.0);
    g = rssringoccs_ThreeVector_Normalize(g);
    a = rssringoccs_Cross_Product(
        rssringoccs_ThreeVector_Rect(0.0, 0.0, 1.0), g);

    a = rssringoccs_ThreeVector_Normalize(a);
    a = rssringoccs_ThreeVector_Scale(0.002, a);
    b = rssringoccs_Cross_Product(g, a);
    b = rssringoccs_ThreeVector_Normalize(b);
    b = rssringoccs_ThreeVector_Scale(0.002, b);

    c = rssringoccs_ThreeVector_Add(a, b);
    c = rssringoccs_ThreeVector_Scale(-256.0, c);
    c = rssringoccs_ThreeVector_Add(c, g);

    fp = fopen("initials_graphic.ppm", "w");
    fprintf(fp, "P6\n%u %u\n255\n", size, size);

    for (y=0; y<size; ++y)
    {
        for(x=0; x<size; ++x)
        {
            p = rssringoccs_ThreeVector_Rect(13.0, 13.0, 13.0);

            for (r=0; r<64; ++r
            )
            {
                factor = (random_unit_interval()-0.5);
                temp1 = rssringoccs_ThreeVector_Scale(factor, a);

                factor = (random_unit_interval()-0.5);
                temp2 = rssringoccs_ThreeVector_Scale(factor, b);

                t = rssringoccs_ThreeVector_Add(temp1, temp2);
                t = rssringoccs_ThreeVector_Scale(99.0, t);

                temp1 = rssringoccs_ThreeVector_Rect(18.0, 20.0, 8.0);
                alpha = rssringoccs_ThreeVector_Add(temp1, t);

                factor = (random_unit_interval() + 0.5*(size-x));
                temp1 = rssringoccs_ThreeVector_Scale(factor, a);

                factor = (random_unit_interval() + 0.5*(size-y));
                temp2 = rssringoccs_ThreeVector_Scale(factor, b);
                beta = rssringoccs_ThreeVector_Add(temp1, temp2);
                beta = rssringoccs_ThreeVector_Add(beta, c);
                beta = rssringoccs_ThreeVector_Scale(16.0, beta);
                temp1 = rssringoccs_ThreeVector_Scale(-1.0, t);
                beta = rssringoccs_ThreeVector_Add(temp1, beta);
                beta = rssringoccs_ThreeVector_Normalize(beta);

                temp1 = Shader(alpha, beta);
                temp2 = rssringoccs_ThreeVector_Scale(3.5, temp1);
                p = rssringoccs_ThreeVector_Add(temp2, p);
            }
            px = rssringoccs_ThreeVector_X(p);
            py = rssringoccs_ThreeVector_Y(p);
            pz = rssringoccs_ThreeVector_Z(p);
            fprintf(fp, "%c%c%c", (int)px, (int)py, (int)pz);
        }
        printf("%u %u\n", y, size);
    }
    return 0;
}