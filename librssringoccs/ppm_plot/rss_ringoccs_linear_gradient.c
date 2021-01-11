#include <rss_ringoccs/include/rss_ringoccs_ppm_plot.h>

void
rssringoccs_RGB_Linear_Gradient(double val, double min, double max, FILE *fp)
{
    double scaled, temp;
    unsigned char red, green, blue;

    if (max < min)
    {
        temp = max;
        max = min;
        min = temp;
    }

    scaled = 255.0*(val - min)/(max-min);

    if (scaled < 64)
    {
        red   = (unsigned char)0;
        green = (unsigned char)(4.0*scaled);
        blue  = (unsigned char)255;
    }
    else if (scaled < 128)
    {
        red   = (unsigned char)0;
        green = (unsigned char)255;
        blue  = (unsigned char)(255 - 4*(scaled - 64));
    }
    else if (scaled < 192)
    {
        red   = (unsigned char)(4.0*(scaled-128.0));
        green = (unsigned char)255;
        blue  = (unsigned char)0;
    }
    else if (scaled < 255)
    {
        red   = (unsigned char)255;
        green = (unsigned char)(255 - 4*(scaled-192));
        blue  = (unsigned char)0;
    }
    else
    {
        red   = (unsigned char)255;
        green = (unsigned char)0;
        blue  = (unsigned char)0;
    }

    /*  Color the current pixel.                                              */
    rssringoccs_Color(red, green, blue, fp);
}
