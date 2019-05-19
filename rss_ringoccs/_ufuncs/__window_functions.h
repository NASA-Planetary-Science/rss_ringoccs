#ifndef RSS_RINGOCCS_WINDOW_FUNCTIONS_H
#define RSS_RINGOCCS_WINDOW_FUNCTIONS_H

#include "__kaiser_bessel.h"

static void __rect(double* wfunc, double w_width, double dx, long nw_pts)
{
        long i;

        for (i=0; i<nw_pts; i++){
            wfunc[i] = 1.0;
        }
}

static void __coss(double* wfunc, double w_width, double dx, long nw_pts)
{
        long i;
        double x;
        dx = ONE_PI * dx / w_width;

        for (i=0; i<nw_pts; i++){
            x = (i-nw_pts) * dx;
            x = cos(x);
            x *= x;
            wfunc[i] = x;
        }
}

static void __kb20(double* wfunc, double w_width, double dx, long nw_pts)
{
    long i;
    double x;

    for (i=0; i<nw_pts; i++){
        x = (i-nw_pts) * dx;
        wfunc[i] = Kaiser_Bessel_Window_2_0(x, w_width);
    }
}

static void __kb25(double* wfunc, double w_width, double dx, long nw_pts)
{
    long i;
    double x;

    for (i=0; i<nw_pts; i++){
        x = (i-nw_pts) * dx;
        wfunc[i] = Kaiser_Bessel_Window_2_5(x, w_width);
    }
}

static void __kb35(double* wfunc, double w_width, double dx, long nw_pts)
{
    long i;
    double x;

    for (i=0; i<nw_pts; i++){
        x = (i-nw_pts) * dx;
        wfunc[i] = Kaiser_Bessel_Window_3_5(x, w_width);
    }
}

static void __kbmd20(double* wfunc, double w_width, double dx, long nw_pts)
{
    long i;
    double x;

    for (i=0; i<nw_pts; i++){
        x = (i-nw_pts) * dx;
        wfunc[i] = Modified_Kaiser_Bessel_Window_2_0(x, w_width);
    }
}

static void __kbmd25(double* wfunc, double w_width, double dx, long nw_pts)
{
    long i;
    double x;

    for (i=0; i<nw_pts; i++){
        x = (i-nw_pts) * dx;
        wfunc[i] = Modified_Kaiser_Bessel_Window_2_5(x, w_width);
    }
}

#endif