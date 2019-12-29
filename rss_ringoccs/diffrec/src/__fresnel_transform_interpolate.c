#include "__diffraction_functions.h"

complex double
Fresnel_Transform_Quadratic_Double(double *x_arr, double *phi_arr,
                                   complex double *T_in, double *w_func,
                                   double kD, double r, double B, double D,
                                   double EPS, long toler, double dx,
                                   double F, long n_pts, long center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double C[3];
    double* psi_n = (double *)malloc(sizeof(double)*n_pts);
    double psi, phi, width;
    double psi_half_mean, psi_full_mean;
    complex double T_out, exp_psi;

    width = x_arr[n_pts-1] - x_arr[0];
    double rcpr_w = 1.0 / width;
    double rcpr_w_sq = rcpr_w * rcpr_w;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = 0.0;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    j = center-(int)((n_pts-1)/2);

     /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){

        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = Newton_Raphson_Fresnel_Psi(kD, r, x_arr[i], phi_arr[i],
                                         phi_arr[i], B, D, EPS, toler);

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi_n[i] = Fresnel_Psi_Double(kD, r, x_arr[i], phi, phi_arr[i], B, D);
    }

    psi_half_mean = (psi_n[(n_pts-1)/4] + psi_n[3*(n_pts-1)/4]) / 2.0;
    psi_full_mean = (psi_n[0] + psi_n[n_pts-1]) / 2;

    free(psi_n);

    C[0] = (16*psi_half_mean-psi_full_mean)*rcpr_w_sq*1.33333333333333333333333;

    for (i = 0; i<n_pts; ++i){

        psi = (r-x_arr[i])*(r-x_arr[i])*C[0];

        exp_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];
        T_out += exp_psi * T_in[j];
        j     += 1;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    T_out *= (0.5 + 0.5*_Complex_I) * dx / F;
    return T_out;
}

complex double
Fresnel_Transform_Quadratic_Norm_Double(double *x_arr, double *phi_arr,
                                        complex double *T_in, double *w_func,
                                        double kD, double r, double B, double D,
                                        double EPS, long toler, double dx,
                                        double F, long n_pts, long center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double C[3];
    double* psi_n = (double *)malloc(sizeof(double)*n_pts);
    double psi, phi, width;
    double psi_half_mean, psi_full_mean;
    complex double T_out, exp_psi, norm;

    width = x_arr[n_pts-1] - x_arr[0];
    double rcpr_w = 1.0 / width;
    double rcpr_w_sq = rcpr_w * rcpr_w;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = 0.0;
    norm  = 0.0;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    j = center-(int)((n_pts-1)/2);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){

        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = Newton_Raphson_Fresnel_Psi(kD, r, x_arr[i], phi_arr[i],
                                         phi_arr[i], B, D, EPS, toler);

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi_n[i] = Fresnel_Psi_Double(kD, r, x_arr[i], phi, phi_arr[i], B, D);
    }

    psi_half_mean = (psi_n[(n_pts-1)/4] + psi_n[3*(n_pts-1)/4]) / 2.0;
    psi_full_mean = (psi_n[0] + psi_n[n_pts-1]) / 2;

    free(psi_n);

    C[0] = (16*psi_half_mean-psi_full_mean)*rcpr_w_sq*1.33333333333333333333333;

    for (i = 0; i<n_pts; ++i){

        psi = (r-x_arr[i])*(r-x_arr[i])*C[0];

        exp_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute the norm using a Riemann sum as well.                     */
        norm   += exp_psi;

        T_out += exp_psi * T_in[j];
        j     += 1;
    }

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number (defined in complex.h). */
    norm = SQRT_2 / cabs(norm);

    /*  Multiply result by the coefficient found in the Fresnel inverse.      *
     *  The 1/F term is omitted, since the F in the norm cancels this.        */
    T_out *= (0.5 + 0.5*_Complex_I) * norm;
    return T_out;
}

complex double
Fresnel_Transform_Cubic_Double(double *x_arr, double *phi_arr,
                               complex double *T_in, double *w_func,
                               double kD, double r, double B, double D,
                               double EPS, long toler, double dx,
                               double F, long n_pts, long center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double C[3];
    double* psi_n = (double *)malloc(sizeof(double)*n_pts);
    double psi, phi, width;
    double psi_half_mean, psi_half_diff;
    double psi_full_mean, psi_full_diff;
    complex double T_out, exp_psi;

    width = x_arr[n_pts-1] - x_arr[0];
    double rcpr_w = 1.0 / width;
    double rcpr_w_sq = rcpr_w * rcpr_w;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = 0.0;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    j = center-(int)((n_pts-1)/2);

     /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){

        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = Newton_Raphson_Fresnel_Psi(kD, r, x_arr[i], phi_arr[i],
                                         phi_arr[i], B, D, EPS, toler);

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi_n[i] = Fresnel_Psi_Double(kD, r, x_arr[i], phi, phi_arr[i], B, D);
    }

    psi_half_mean = (psi_n[(n_pts-1)/4] + psi_n[3*(n_pts-1)/4]) / 2.0;
    psi_full_mean = (psi_n[0] + psi_n[n_pts-1]) / 2;
    psi_half_diff = psi_n[(n_pts-1)/4] - psi_n[3*(n_pts-1)/4];
    psi_full_diff = psi_n[0] - psi_n[n_pts-1];

    free(psi_n);

    C[0] = (16*psi_half_mean-psi_full_mean)*rcpr_w_sq*1.33333333333333333333333;
    C[1] = (psi_full_diff-2.0*psi_half_diff)*rcpr_w_sq*rcpr_w*5.333333333333333;
    for (i = 0; i<n_pts; ++i){

        psi = C[1];
        psi = psi*(r-x_arr[i])+C[0];
        psi = psi*(r-x_arr[i])*(r-x_arr[i]);

        exp_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];
        T_out += exp_psi * T_in[j];
        j     += 1;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    T_out *= (0.5 + 0.5*_Complex_I) * dx / F;
    return T_out;
}

complex double
Fresnel_Transform_Cubic_Norm_Double(double *x_arr, double *phi_arr,
                                    complex double *T_in, double *w_func,
                                    double kD, double r, double B, double D,
                                    double EPS, long toler, double dx,
                                    double F, long n_pts, long center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double C[3];
    double* psi_n = (double *)malloc(sizeof(double)*n_pts);
    double psi, phi, width;
    double psi_half_mean, psi_half_diff;
    double psi_full_mean, psi_full_diff;
    complex double T_out, exp_psi, norm;

    width = x_arr[n_pts-1] - x_arr[0];
    double rcpr_w = 1.0 / width;
    double rcpr_w_sq = rcpr_w * rcpr_w;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = 0.0;
    norm  = 0.0;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    j = center-(int)((n_pts-1)/2);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){

        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = Newton_Raphson_Fresnel_Psi(kD, r, x_arr[i], phi_arr[i],
                                         phi_arr[i], B, D, EPS, toler);

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi_n[i] = Fresnel_Psi_Double(kD, r, x_arr[i], phi, phi_arr[i], B, D);
    }

    psi_half_mean = (psi_n[(n_pts-1)/4] + psi_n[3*(n_pts-1)/4]) / 2.0;
    psi_full_mean = (psi_n[0] + psi_n[n_pts-1]) / 2;
    psi_half_diff = psi_n[(n_pts-1)/4] - psi_n[3*(n_pts-1)/4];
    psi_full_diff = psi_n[0] - psi_n[n_pts-1];

    free(psi_n);

    C[0] = (16*psi_half_mean-psi_full_mean)*rcpr_w_sq*1.33333333333333333333333;
    C[1] = (psi_full_diff-2.0*psi_half_diff)*rcpr_w_sq*rcpr_w*5.333333333333333;

    for (i = 0; i<n_pts; ++i){

        psi = C[1];
        psi = psi*(r-x_arr[i])+C[0];
        psi = psi*(r-x_arr[i])*(r-x_arr[i]);

        exp_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute the norm using a Riemann sum as well.                     */
        norm   += exp_psi;

        T_out += exp_psi * T_in[j];
        j     += 1;
    }

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number (defined in complex.h). */
    norm = SQRT_2 / cabs(norm);

    /*  Multiply result by the coefficient found in the Fresnel inverse.      *
     *  The 1/F term is omitted, since the F in the norm cancels this.        */
    T_out *= (0.5 + 0.5*_Complex_I) * norm;
    return T_out;
}

complex double
Fresnel_Transform_Quartic_Double(double *x_arr, double *phi_arr,
                                 complex double *T_in, double *w_func,
                                 double kD, double r, double B, double D,
                                 double EPS, long toler, double dx,
                                 double F, long n_pts, long center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double C[3];
    double* psi_n = (double *)malloc(sizeof(double)*n_pts);
    double psi, phi, width;
    double psi_half_mean, psi_half_diff;
    double psi_full_mean, psi_full_diff;
    complex double T_out, exp_psi;

    width = x_arr[n_pts-1] - x_arr[0];
    double rcpr_w = 1.0 / width;
    double rcpr_w_sq = rcpr_w * rcpr_w;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = 0.0;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    j = center-(int)((n_pts-1)/2);

     /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){

        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = Newton_Raphson_Fresnel_Psi(kD, r, x_arr[i], phi_arr[i],
                                         phi_arr[i], B, D, EPS, toler);

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi_n[i] = Fresnel_Psi_Double(kD, r, x_arr[i], phi, phi_arr[i], B, D);
    }

    psi_half_mean = (psi_n[(n_pts-1)/4] + psi_n[3*(n_pts-1)/4]) / 2.0;
    psi_full_mean = (psi_n[0] + psi_n[n_pts-1]) / 2;
    psi_half_diff = psi_n[(n_pts-1)/4] - psi_n[3*(n_pts-1)/4];
    psi_full_diff = psi_n[0] - psi_n[n_pts-1];

    free(psi_n);

    C[0] = (16*psi_half_mean-psi_full_mean)*rcpr_w_sq*1.33333333333333333333333;
    C[1] = (psi_full_diff-2.0*psi_half_diff)*rcpr_w_sq*rcpr_w*5.333333333333333;
    C[2] = (psi_full_mean-4.0*psi_half_mean)*rcpr_w_sq*rcpr_w_sq*21.33333333333;

    for (i = 0; i<n_pts; ++i){

        psi = C[2];
        psi = psi*(r-x_arr[i])+C[1];
        psi = psi*(r-x_arr[i])+C[0];
        psi = psi*(r-x_arr[i])*(r-x_arr[i]);

        exp_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];
        T_out += exp_psi * T_in[j];
        j     += 1;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    T_out *= (0.5 + 0.5*_Complex_I) * dx / F;
    return T_out;
}

complex double
Fresnel_Transform_Quartic_Norm_Double(double *x_arr, double *phi_arr,
                                      complex double *T_in, double *w_func,
                                      double kD, double r, double B, double D,
                                      double EPS, long toler, double dx,
                                      double F, long n_pts, long center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    long i, j;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double C[3];
    double* psi_n = (double *)malloc(sizeof(double)*n_pts);
    double psi, phi, width;
    double psi_half_mean, psi_half_diff;
    double psi_full_mean, psi_full_diff;
    complex double T_out, exp_psi, norm;

    width = x_arr[n_pts-1] - x_arr[0];
    double rcpr_w = 1.0 / width;
    double rcpr_w_sq = rcpr_w * rcpr_w;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = 0.0;
    norm  = 0.0;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    j = center-(int)((n_pts-1)/2);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (i = 0; i<n_pts; ++i){

        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = Newton_Raphson_Fresnel_Psi(kD, r, x_arr[i], phi_arr[i],
                                         phi_arr[i], B, D, EPS, toler);

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        psi_n[i] = Fresnel_Psi_Double(kD, r, x_arr[i], phi, phi_arr[i], B, D);
    }

    psi_half_mean = (psi_n[(n_pts-1)/4] + psi_n[3*(n_pts-1)/4]) / 2.0;
    psi_full_mean = (psi_n[0] + psi_n[n_pts-1]) / 2;
    psi_half_diff = psi_n[(n_pts-1)/4] - psi_n[3*(n_pts-1)/4];
    psi_full_diff = psi_n[0] - psi_n[n_pts-1];

    free(psi_n);

    C[0] = (16*psi_half_mean-psi_full_mean)*rcpr_w_sq*1.33333333333333333333333;
    C[1] = (psi_full_diff-2.0*psi_half_diff)*rcpr_w_sq*rcpr_w*5.333333333333333;
    C[2] = (psi_full_mean-4.0*psi_half_mean)*rcpr_w_sq*rcpr_w_sq*21.33333333333;

    for (i = 0; i<n_pts; ++i){

        psi = C[2];
        psi = psi*(r-x_arr[i])+C[1];
        psi = psi*(r-x_arr[i])+C[0];
        psi = psi*(r-x_arr[i])*(r-x_arr[i]);

        exp_psi = (cos(psi) - _Complex_I*sin(psi))*w_func[i];

        /*  Compute the norm using a Riemann sum as well.                     */
        norm   += exp_psi;

        T_out += exp_psi * T_in[j];
        j     += 1;
    }

    /*  The integral in the numerator of norm evaluates to F sqrt(2). Use     *
     *  this in the calculation of the normalization. The cabs function       *
     *  computes the absolute value of complex number (defined in complex.h). */
    norm = SQRT_2 / cabs(norm);

    /*  Multiply result by the coefficient found in the Fresnel inverse.      *
     *  The 1/F term is omitted, since the F in the norm cancels this.        */
    T_out *= (0.5 + 0.5*_Complex_I) * norm;
    return T_out;
}