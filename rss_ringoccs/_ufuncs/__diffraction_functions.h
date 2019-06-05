#ifndef RSS_RINGOCCS_DIFFRACTION_FUNCTIONS_H
#define RSS_RINGOCCS_DIFFRACTION_FUNCTIONS_H

#include <math.h>
#include <complex.h>

/*  Various coefficients and constants defined here.      */
#include "__math_constants.h"

static void get_arr(double* x_arr, double dx, long nw_pts)
{
    /***************************************************************************
     *  Function:                                                              *
     *      get_arr                                                            *
     *  Purpose:                                                               *
     *      This computes the array pi/2 x^2, where x range from rho-w/2 to    *
     *      rho+w/2. Do the the symmetry involved in the computation of the    *
     *      Fresnel transform, on the points from rho-w/2 <= x < 0 are         *
     *      computed. The midpoint is always 0, and the right side is equal to *
     *      the left side. Hence, given a window with 2n+1 points, the array   *
     *      passed to get_arr will have n points. This cuts the number of      *
     *      computations in half while computing the same result.              *
     **************************************************************************/
    long i;
    for (i=0; i<nw_pts; ++i){
        x_arr[i] = (i-nw_pts)*dx;
    }
}

complex double _fresnel_transform(double* x_arr, char* T_in, double* w_func,
                                  double F, double dx, long n_pts,
                                  npy_intp T_in_steps)
{
    long i, j;
    double x, F2;
    complex double T_out, exp_negative_ix;

    T_out = 0.0;

    j = -n_pts;

    /*  Division is more expensive than multiplication, so store the
        reciprical of F as a variable and compute with that.          */
    F = 1.0/F;
    F2 = F*F;

    for (i = 0; i<n_pts; ++i){
        x = x_arr[i]*F2;
        exp_negative_ix = (cos(x) - _Complex_I*sin(x)) * w_func[i];
        T_out += exp_negative_ix * (*(complex double *)(T_in + j*T_in_steps) +
                                    *(complex double *)(T_in - j*T_in_steps));
        j += 1;
    }
        T_out += *(complex double *)T_in;
        T_out *= (0.5+0.5*_Complex_I)*dx*F;
    return T_out;
}

complex double _fresnel_transform_norm(double* x_arr, char* T_in,
                                       double* w_func, double F, double dx,
                                       long n_pts, npy_intp T_in_steps)
{
    long i, j;
    double x, F2;
    complex double T_out, exp_negative_ix, norm;

    T_out  = 0.0;
    norm   = 0.0;

    j = -n_pts;

    /*  Division is more expensive than multiplication, so store the
        reciprical of F as a variable and compute with that.          */
    F = 1.0/F;
    F2 = F*F;

    for (i = 0; i<n_pts; ++i){
        x = x_arr[i]*F2;
        exp_negative_ix = (cos(x) - _Complex_I*sin(x)) * w_func[i];
        norm  += 2.0*exp_negative_ix;
        T_out += exp_negative_ix * (*(complex double *)(T_in + j*T_in_steps) +
                                    *(complex double *)(T_in - j*T_in_steps));
        j += 1;
    }
        norm  += 1.0;
        norm   = SQRT_2 / (F*cabs(norm));
        T_out += *(complex double *)T_in;
        T_out *= (0.5+0.5*_Complex_I)*F*norm;
    return T_out;
}

complex double _fresnel_quartic(double* x_arr, char* T_in, double* w_func,
                                double rcpr_D, double A_0, double A_1,
                                double A_2, double dx, double rcpr_F, double kd,
                                long n_pts, npy_intp T_in_steps)
{
    /***************************************************************************
     *  Compute the Fresnel Transform:                                         *
     *                                                                         *
     *              W/2                                                        *
     *               -                                                         *
     *              | |                                                        *
     *   T(rho) =   |   T_hat(rho_0)w(rho-rho_0)exp(-i psi) drho_0             *
     *            | |                                                          *
     *             -                                                           *
     *            -W/2                                                         *
     *                                                                         *
     *  Where W is the window width, psi is the Fresnel Kernel, T_hat is the   *
     *  diffrected complex transmittance, and w is the window function.        *
     **************************************************************************/
    long i, j;
    double x, x2, psi, psi_even, psi_odd;
    complex double T_out, exp_negative_psi, exp_positive_psi;

    /*  Initialize T_out to zero.  */
    T_out = 0.0;
    j = -n_pts;

    for (i = 0; i<n_pts; ++i){
        x   = x_arr[i]*rcpr_D;
        x2  = x*x;

        psi_even    = kd*x2*(A_0+x2*A_2);
        psi_odd     = kd*x2*x*A_1;

        /*  Compute exp(-ipsi) using Euler's Formula.  */
        psi = psi_even - psi_odd;
        exp_negative_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        psi = psi_even + psi_odd;
        exp_positive_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Approximate the integral with a Riemann Sum.  */
        T_out += exp_negative_psi * *(complex double *)(T_in + j*T_in_steps);
        T_out += exp_positive_psi * *(complex double *)(T_in - j*T_in_steps);
        j += 1;
    }
        T_out += *(complex double *)T_in;
        T_out *= (0.5 + 0.5*_Complex_I) * dx * rcpr_F;
    return T_out;
}

#endif