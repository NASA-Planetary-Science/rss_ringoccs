/******************************************************************************
 *                          Diffraction Functions                             *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      This file contains routines to compute the Fresnel transform for      *
 *      diffraction reconstruction or forward modeling.                       *
 *          Fresnel Quadratic Approximations:                                 *
 *              Classic quadratic approximation used in Fourier optics.       *
 *          Legendre Expansion:                                               *
 *             Approximate the Fresnel kernel using Legendre polynomials.     *
 *          Newton-Raphson Method:                                            *
 *              Performs the Fresnel Transform by computing the stationary    *
 *              value of the Fresnel Kernel using the Newton-Raphson method.  *
 ******************************************************************************
 *  The Inverse Fresnel Transform:                                            *
 *                                                                            *
 *                infinity                                                    *
 *                     -                                                      *
 *                    | |                                                     *
 *         T(rho) =   |   T_hat(r_0)w(r-r_0)exp(-i psi(r,r_0)) dr_0           *
 *                  | |                                                       *
 *                   -                                                        *
 *               -infinity                                                    *
 *                                                                            *
 *      Where T_hat is the diffracted data, w is the window function, r is    *
 *      the ring intercept point, and r_0 is a dummy variable of integration. *
 *      psi is the Fresnel Kernel, and exp is the exponential function.       *
 ******************************************************************************
 *  The Normalization Scheme:                                                 *
 *      As the resolution get's too high, say 10 km or larger, the window     *
 *      width quickly shrinks to zero and the integral will be approximately  *
 *      zero. To account for this, the option to normalize the integral by    *
 *      the window width is offered. The normalization is defined as follows: *
 *                                                                            *
 *                    |     _ +infinity           |                           *
 *                    |    | |                    |                           *
 *                    |    |    exp(-i psi(x)) dx |                           *
 *                    |  | |                      |                           *
 *                    |   -  -infinity            |                           *
 *          Norm =  __________________________________                        *
 *                  |    -  +W/2                    |                         *
 *                  |   | |                         |                         *
 *                  |   |    w(x) exp(-i psi(x)) dx |                         *
 *                  | | |                           |                         *
 *                  |  -   -W/2                     |                         *
 *                                                                            *
 *      This has the effect of making the free-space regions, which are       *
 *      regions that are not affected by diffraction, evaluate to one,        *
 *      regardless of what (positive) resolution was chosen.                  *
 ******************************************************************************
 *                              DEFINED FUNCTIONS                             *
 ******************************************************************************
 *  reset_window:                                                             *
 *      Void function that takes in a pointer to a double array and creates   *
 *      an array of values for the ring radius within half a window width of  *
 *      the ring intercept point. Do to symmetry, only the values to the left *
 *      of the ring intercept point are computed. This also computes the      *
 *      window function as a function of this array.                          *
 ******************************************************************************
 *  Fresnel_Transform_Double / Fresnel_Transform_Norm_Double                  *
 *      Computes the Fresnel transform using a quadratic approximation to the *
 *      Fresnel kernel, with or without normalization.                        *
 ******************************************************************************
 *  Fresnel_Legendre_Double / Fresnel_Legendre_Norm_Double                    *
 *      Computes the Fresnel Transform using Legendre Polynomials to          *
 *      approximate the Fresnel kernel to various powers (4, 6, or 8).        *
 ******************************************************************************
 *  Fresnel_Transform_Newton_Double / Fresnel_Transform_Newton_Norm_Double    *
 *      Computes the Fresnel inverse transform using Newton-Raphson to        *
 *      compute the stationary value of the Fresnel kernel.                   *
 ******************************************************************************
 *                            A FRIENDLY WARNING                              *
 ******************************************************************************
 *  This code uses complex numbers throughout, and is compatible with the C99 *
 *  standard. To use this code, make sure your compiler supports C99 or more  *
 *  recent standards of the C Programming Language.                           *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       June 21, 2019                                                 *
 ******************************************************************************/
#include "__diffraction_functions.h"

/******************************************************************************
 *  Function:                                                                 *
 *      reset_window                                                          *
 *  Purpose:                                                                  *
 *      Compute an array of points from -width to zero, equally spaced by dx. *
 *      This acts as the independent variable for later use.                  *
 *  Arguments:                                                                *
 *      x_arr (double *):                                                     *
 *          Defined as rho-rho0, where rho0 is the ring radius of the point   *
 *          being reconstructed, and rho is the dummy variable of integration *
 *          which varies from rho0-W/2 to rho0, W being the window width.     *
 *      w_func (double *):                                                    *
 *          The window/tapering function, as a function of x_arr.             *
 *      dx (double):                                                          *
 *          The sample spacing, equivalent to x_arr[1] - x_arr[0].            *
 *      width (double):                                                       *
 *          The width of the window function.                                 *
 *      nw_pts (long):                                                        *
 *          Half the number of points in the window width. The symmetry of    *
 *          the quadratic approximation allows one to perform the inversion   *
 *          with only half of the window. This saves a lot of computation.    *
 *      fw  double (*)(double, double):                                       *
 *          Function pointer to the window function.                          *
 *  Notes:                                                                    *
 *      1.) This is a void function that takes in pointers as arguments. The  *
 *          values of the pointers are changed within this function and there *
 *          is no need to return anything. Hence, no return statement.        *
 ******************************************************************************/
static void reset_window(double *x_arr, double *w_func, double dx, double width,
                         long nw_pts, double (*fw)(double, double))
{
    /*  Create a variable for indexing.                                       */
    long i;

    /* Loop over i, computing the window function and the x_arr variable.     */
    for(i=0; i<nw_pts; ++i){
        x_arr[i] = (i-nw_pts)*dx;
        w_func[i] = fw(x_arr[i], width);
    }
}

static int check_dlp_data(DLPObj *dlp)
{
    /*  If any of these pointers are void, return 0. Else, pass.              */
    if ((!dlp->T_in) || (!dlp->rho_km_vals)  || (!dlp->F_km_vals)
                     || (!dlp->phi_rad_vals) || (!dlp->kd_vals)
                     || (!dlp->B_rad_vals)   || (!dlp->D_km_vals)
                     || (!dlp->w_km_vals)    || (!dlp->T_out)) return 0;
    else return 1;
}

void DiffractionCorrectionFresnel(DLPObj *dlp)
{
    /*  If everything executes smoothly, status should remain at zero.        */
    dlp->status = 0;

    /*  Check that the pointers to the data are not NULL.                     */
    if (!(check_dlp_data(dlp)))
    {
        /*  One of the variables has null data, return to calling function.   */
        dlp->status = 1;
        return;
    }

    /*  i and j used for indexing, nw_pts is number of points in window.     */
    long i, j, nw_pts, center;

    /*  w_init is window width (km), dx and two_dx are sample spacing (km).  */
    double w_init, dx, two_dx;

    /*  Declare function pointers. fw is window function, FresT is Fresnel   *
     *  transform. Options for FresT are with or without normalization.      */
    double (*fw)(double, double);
    complex double (*FresT)(double*, complex double*, double*, double,
                            double, long, long);

    /*  Cast the selected window type to the fw pointer.                     */
    if      (dlp->wtype == 0) fw = &Rect_Window_Double;
    else if (dlp->wtype == 1) fw = &Coss_Window_Double;
    else if (dlp->wtype == 2) fw = &Kaiser_Bessel_2_0_Double;
    else if (dlp->wtype == 3) fw = &Kaiser_Bessel_2_5_Double;
    else if (dlp->wtype == 4) fw = &Kaiser_Bessel_3_5_Double;
    else if (dlp->wtype == 5) fw = &Modified_Kaiser_Bessel_2_0_Double;
    else if (dlp->wtype == 6) fw = &Modified_Kaiser_Bessel_2_5_Double;
    else                      fw = &Modified_Kaiser_Bessel_3_5_Double;

    /*  Cast FresT to the appropriate function.                               */
    if (dlp->use_norm) FresT = &Fresnel_Transform_Norm_Double;
    else               FresT = &Fresnel_Transform_Double;

    /* Move the pointers to the correct starting point. Compute window width.*/
    center = dlp->start;

    /*  Compute some extra necessary variables.                               */
    w_init  = dlp->w_km_vals[center];
    dx      = dlp->rho_km_vals[center+1] - dlp->rho_km_vals[center];
    two_dx  = 2.0*dx;
    nw_pts  = ((long)(w_init / two_dx))+1;

    if (center - nw_pts < 0)
    {
        /*  Invalid index for the for the first point. That is, the window    *
         *  of integration required for the first point goes beyond the       *
         *  available data. Return to calling function.                       */
        dlp->n_used = nw_pts;
        dlp->status = 2;
        return;
    }

    /*  Reserve some memory for two arrays, the ring radius and the window    *
     *  function. This will need to be reallocated later if the window width  *
     *  changes by more than two_dx.                                          */
    double* x_arr  = (double *)malloc(sizeof(double)*nw_pts);
    double* w_func = (double *)malloc(sizeof(double)*nw_pts);


    /*  Check that malloc was successfull then pass the x_arr array           *
     *  (ring radius) to the void function reset_window. This alters x_arr so *
     *  that it's values range from -W/2 to zero, W begin the window width.   */
    if (!(x_arr) || !(w_func))
    {
        /*  Malloc failed, return to calling function.                        */
        dlp->status = 3;
        return;
    }
    else reset_window(x_arr, w_func, dx, w_init, nw_pts, fw);

    /* Compute Window Functions, and compute pi/2 * x^2                       */
    for(j=0; j<nw_pts; ++j)
    {
        /*  The independent variable is pi/2 * ((rho-rho0)/F)^2. Compute      *
         *  part of this. The 1/F^2 part is introduced later.                 */
        x_arr[j] *= PI_BY_TWO*x_arr[j];

        /*  If forward transform is selected, negate x_arr.                   */
        x_arr[j] *= (dlp->use_fwd == 0) - (dlp->use_fwd == 1);
    }

    /*  Compute the Fresnel transform across the input data.                  */
    for (i=0; i<=dlp->n_used; ++i)
    {
        /*  If the window width has deviated more the 2*dx, reset variables.  */
        if (fabs(w_init - dlp->w_km_vals[center]) >= two_dx)
        {
            /* Reset w_init and recompute window function.                    */
            w_init = dlp->w_km_vals[center];
            nw_pts = ((long)(w_init / two_dx))+1;

            /*  Reallocate memory, since the sizes of the arrays changed.     */
            w_func = (double *)realloc(w_func, sizeof(double)*nw_pts);
            x_arr  = (double *)realloc(x_arr, sizeof(double)*nw_pts);

            /*  Reset the x_arr array to range between -W/2 and zero.         */
            reset_window(x_arr, w_func, dx, w_init, nw_pts, fw);

            /* Compute Window Functions, and compute pi/2 * x^2               */
            for(j=0; j<nw_pts; ++j)
            {
                x_arr[j] *= PI_BY_TWO*x_arr[j];

                /*  Again, if forward calculation is set, negate x_arr.       */
                x_arr[j] *= (dlp->use_fwd == 0) - (dlp->use_fwd == 1);
            }
        }

        /*  Compute the Fresnel Transform about the current point.            */
        dlp->T_out[i] = FresT(x_arr, dlp->T_in, w_func,
                              dlp->F_km_vals[center], dx, nw_pts, center);

        /*  Move the pointers to the next point.                              */
        center += 1;
    }

    /*  Free the variables allocated by malloc.                               */
    free(x_arr);
    free(w_func);
}

void DiffractionCorrectionLegendre(DLPObj *dlp)
{
    /*  If everything executes smoothly, status should remain at zero.        */
    dlp->status = 0;

    /*  Check that the pointers to the data are not NULL.                     */
    if (!(check_dlp_data(dlp)))
    {
        /*  One of the variables has null data, return to calling function.   */
        dlp->status = 1;
        return;
    }

    /*  i is used for indexing, nw_pts is the number of points in the window. */
    long i, nw_pts, center;

    /*  IsEven is a boolean for determining the parity of the polynomial.     */
    bool IsEven;

    /*  Variable for the number of Legendre coefficients to be computed.      */
    unsigned char poly_order;

    /*  Various other variables needed throughout.                            */
    double w_init, dx, two_dx, cosb, sinp, cosp, Legendre_Coeff;

    /*  Create function pointers for window function and Fresnel transform.   */
    double (*fw)(double, double);
    complex double (*FresT)(double*, complex double*, double*, double, double *,
                            double, double, double, long, unsigned char, long);

    /*  Cast the selected window type to the fw pointer.                      */
    if      (dlp->wtype == 0) fw = &Rect_Window_Double;
    else if (dlp->wtype == 1) fw = &Coss_Window_Double;
    else if (dlp->wtype == 2) fw = &Kaiser_Bessel_2_0_Double;
    else if (dlp->wtype == 3) fw = &Kaiser_Bessel_2_5_Double;
    else if (dlp->wtype == 4) fw = &Kaiser_Bessel_3_5_Double;
    else if (dlp->wtype == 5) fw = &Modified_Kaiser_Bessel_2_0_Double;
    else if (dlp->wtype == 6) fw = &Modified_Kaiser_Bessel_2_5_Double;
    else                      fw = &Modified_Kaiser_Bessel_3_5_Double;

    /*  Set the IsEven boolean to the appropriate value. Since the linear and *
     *  constant term are zero, even polynomials will have an odd number of   *
     *  terms. For example, the quartic expansion has the quadratic, cubic,   *
     *  quartic terms, and hence needs three coefficients. Thus, if the       *
     *  dlp->order variable is an odd number, this corresponds to an even     *
     *  polynomial and vice versa. Set IsEven accordingly.                    */
    if (2*(dlp->order/2) == dlp->order) IsEven = False;
    else IsEven = True;

    /*  Select the appropriate Fresnel transform and set poly_order.          */
    if (IsEven)
    {
        poly_order = dlp-> order;

        /*  Select the even transform, with or without normalization.         */
        if (dlp->use_norm) FresT = &Fresnel_Transform_Legendre_Norm_Even_Double;
        else               FresT = &Fresnel_Transform_Legendre_Even_Double;
    }
    else 
    {
        poly_order = dlp->order + 1;

        /*  Choose the odd transform with or without normalization.           */
        if (dlp->use_norm) FresT = &Fresnel_Transform_Legendre_Norm_Odd_Double;
        else               FresT = &Fresnel_Transform_Legendre_Odd_Double;
    }

    /* Compute first window width and window function. */
    center = dlp->start;

    /*  If forward tranform is set, negate the kd_vals variable. This will    *
     *  the equivalent effect of computing the forward calculation later.     */
    if (dlp->use_fwd)
    {
        /*  Loop over all of kd_vals and negate the value.                    */
        for (i=0; i <= dlp->n_used; ++i)
        {
            dlp->kd_vals[i] *= -1.0;
        }
    }

    /*  Compute more necessary data.                                          */
    w_init  = dlp->w_km_vals[center];
    dx      = dlp->rho_km_vals[center+1] - dlp->rho_km_vals[center];
    two_dx  = 2.0*dx;
    nw_pts  = (long)(w_init / two_dx)+1;

    if (center - nw_pts < 0)
    {
        /*  Invalid index for the for the first point. That is, the window    *
         *  of integration required for the first point goes beyond the       *
         *  available data. Return to calling function.                       */
        dlp->n_used = nw_pts;
        dlp->status = 2;
        return;
    }

    /*  Allocate memory for the independent variable and window function.     */
    double *x_arr          = (double *)malloc(sizeof(double)*nw_pts);
    double *w_func         = (double *)malloc(sizeof(double)*nw_pts);

    /*  Also for the two Legendre polynomials.                                */
    double *legendre_p     = (double *)malloc(sizeof(double)*(poly_order+1));
    double *alt_legendre_p = (double *)malloc(sizeof(double)*poly_order);

    /*  And finally for the coefficients of psi.                              */
    double *fresnel_ker_coeffs = (double *)malloc(sizeof(double)*poly_order);

    /*  Check that malloc was successfull then pass the x_arr array           *
     *  (ring radius) to the void function reset_window. This alters x_arr so *
     *  that it's values range from -W/2 to zero, W begin the window width.   */
    if (!(x_arr)    ||    !(w_func)            ||    !(legendre_p)
                    ||    !(alt_legendre_p)    ||    !(fresnel_ker_coeffs))
    {
        /*  Malloc failed, return to calling function.                        */
        dlp->status = 3;
        return;
    }
    else reset_window(x_arr, w_func, dx, w_init, nw_pts, fw);

    for (i = 0; i <= dlp->n_used; ++i)
    {
        /*  Compute some geometric information, and the scaling coefficient   *
         *  for the Legendre polynomial expansion.                            */
        cosb            = cos(dlp->B_rad_vals[center]);
        cosp            = cos(dlp->phi_rad_vals[center]);
        sinp            = sin(dlp->phi_rad_vals[center]);
        Legendre_Coeff  = cosb*sinp;
        Legendre_Coeff *= Legendre_Coeff;
        Legendre_Coeff  = 0.5*Legendre_Coeff/(1.0-Legendre_Coeff);

        /* Compute Legendre Polynomials,                                      */
        Legendre_Polynomials(legendre_p, cosb*cosp, poly_order+1);
        Alt_Legendre_Polynomials(alt_legendre_p, legendre_p, poly_order);

        /*  Compute the coefficients using Cauchy Products. First compute     *
         *  the bottom triangle of the square in the product.                 */
        Fresnel_Kernel_Coefficients(fresnel_ker_coeffs, legendre_p,
                                    alt_legendre_p, Legendre_Coeff, poly_order);

        /*  If the window width changes significantly, recompute w_func.      */
        if (fabs(w_init - dlp->w_km_vals[center]) >= two_dx)
        {
            /* Reset w_init and recompute window function.                    */
            w_init = dlp->w_km_vals[center];
            nw_pts = ((long)(w_init / two_dx))+1;

            /*  Reallocate x_arr and w_func since the sizes changed.          */
            x_arr  = (double *)realloc(x_arr, sizeof(double)*nw_pts);
            w_func = (double *)realloc(w_func, sizeof(double)*nw_pts);

            /*  Recompute x_arr and w_func for the new sizes.                 */
            reset_window(x_arr, w_func, dx, w_init, nw_pts, fw);
        }

        /*  Compute the fresnel tranform about the current point.             */
        dlp->T_out[i] = FresT(x_arr, dlp->T_in, w_func, dlp->D_km_vals[center],
                              fresnel_ker_coeffs, dx, dlp->F_km_vals[center],
                              dlp->kd_vals[center], nw_pts, dlp->order, center);

        /*  Increment T_in pointer using pointer arithmetic.                  */
        center += 1;
    }

    /*  Free all variables allocated by malloc.                               */
    free(x_arr);
    free(w_func);
    free(legendre_p);
    free(alt_legendre_p);
    free(fresnel_ker_coeffs);
}

void DiffractionCorrectionNewton(DLPObj *dlp)
{
    /*  If everything executes smoothly, status should remain at zero.        */
    dlp->status = 0;

    /*  Check that the pointers to the data are not NULL.                     */
    if (!(check_dlp_data(dlp)))
    {
        /*  One of the variables has null data, return to calling function.   */
        dlp->status = 1;
        return;
    }

    /*  Variables for indexing. nw_pts is the number of points in the window. */
    long i, j, nw_pts, center;

    /*  Toler is the number of iterations allowed in Newton-Raphson.          */
    long toler;

    /*  Some variables needed for reconstruction.                             */
    double w_init, dx, two_dx;

    /*  EPS is the maximum allowed error in the Newton-Raphson scheme.        */
    double EPS;

    /*  Set toler to 5 and EPS to e-4, reasonable for almost all cases.       */
    toler = 5;
    EPS = 1.E-4;

    /*  Function pointers for the window function and the Fresnel transform.  */
    double (*fw)(double, double);
    complex double (*FresT)(double *, double *, complex double *, double *,
                            double, double, double, double, double, long,
                            double, double, long, long);

    /*  Cast the selected window type to the fw pointer.                      */
    if      (dlp->wtype == 0) fw = &Rect_Window_Double;
    else if (dlp->wtype == 1) fw = &Coss_Window_Double;
    else if (dlp->wtype == 2) fw = &Kaiser_Bessel_2_0_Double;
    else if (dlp->wtype == 3) fw = &Kaiser_Bessel_2_5_Double;
    else if (dlp->wtype == 4) fw = &Kaiser_Bessel_3_5_Double;
    else if (dlp->wtype == 5) fw = &Modified_Kaiser_Bessel_2_0_Double;
    else if (dlp->wtype == 6) fw = &Modified_Kaiser_Bessel_2_5_Double;
    else                      fw = &Modified_Kaiser_Bessel_3_5_Double;

    /*  Select the correct Fresnel transformation.                            */
    if (dlp->use_norm)
    {
        if (dlp->interp == 0)
        {
            FresT = &Fresnel_Transform_Newton_Norm_Double;
        }
        else if (dlp->interp == 2)
        {
            FresT = &Fresnel_Transform_Quadratic_Norm_Double;
        }
        else if (dlp->interp == 3)
        {
            FresT = &Fresnel_Transform_Cubic_Norm_Double;
        }
        else if (dlp->interp == 4)
        {
            FresT = &Fresnel_Transform_Quartic_Norm_Double;
        }
        else
        {
            /*  Illegal input for interp, return to calling function.         */
            dlp->status = 4;
            return;
        }
    }
    else
    {
        if (dlp->interp == 0)      FresT = &Fresnel_Transform_Newton_Double;
        else if (dlp->interp == 2) FresT = &Fresnel_Transform_Quadratic_Double;
        else if (dlp->interp == 3) FresT = &Fresnel_Transform_Cubic_Double;
        else if (dlp->interp == 4) FresT = &Fresnel_Transform_Quartic_Double;
        else
        {
            /*  Illegal input for interp, return to calling function.         */
            dlp->status = 4;
            return;
        }
    }

    /* Compute first window width and window function. */
    center = dlp->start;

    /*  If forward tranform is set, negate the kd_vals variable. This will    *
     *  the equivalent effect of computing the forward calculation later.     */
    if (dlp->use_fwd)
    {
        /*  Loop over all of kd_vals and negate the value.                    */
        for (i=0; i <= dlp->n_used; ++i)
        {
            dlp->kd_vals[i] *= -1.0;
        }
    }

    /*  Compute some more variables.                                          */
    w_init  = dlp->w_km_vals[center];
    dx      = dlp->rho_km_vals[center+1] - dlp->rho_km_vals[center];
    two_dx  = 2.0*dx;
    nw_pts  = 2*((long)(w_init / (2.0 * dx)))+1;

    if (center - ((nw_pts-1)/2) < 0)
    {
        /*  Invalid index for the for the first point. That is, the window    *
         *  of integration required for the first point goes beyond the       *
         *  available data. Return to calling function.                       */
        dlp->n_used = (nw_pts-1)/2;
        dlp->status = 2;
        return;
    }

    /*  Allocate memory for these required variables.                         */
    double *x_arr   = (double *)malloc(sizeof(double) * nw_pts);
    double *phi_arr = (double *)malloc(sizeof(double) * nw_pts);
    double *w_func  = (double *)malloc(sizeof(double) * nw_pts);

    /*  Check that malloc was successfull.                                    */
    if (!(x_arr) || !(w_func) || !(phi_arr))
    {
        /*  Malloc failed, return to calling function.                        */
        dlp->status = 3;
        return;
    }

    /*  Compute the rho and phi variables, and the window function.           */
    for (j=0; j<nw_pts; ++j)
    {
        x_arr[j]   = dlp->rho_km_vals[center+j-(nw_pts-1)/2];
        phi_arr[j] = dlp->phi_rad_vals[center+j-(nw_pts-1)/2];
        w_func[j]  = fw(x_arr[j] - dlp->rho_km_vals[center], w_init);
    }

    /*  Run diffraction correction point by point.                            */
    for (i=0; i<=dlp->n_used; ++i)
    {
        /*  If the window width changes significantly, recompute w_func.      */
        if (fabs(w_init - dlp->w_km_vals[center]) >= two_dx)
        {
            /* Reset w_init and recompute window function.                    */
            w_init  = dlp->w_km_vals[center];
            nw_pts  = 2*((int)(w_init / (2.0 * dx)))+1;

            /*  Reallocate memory since the sizes have changed.               */
            w_func  = (double *)realloc(w_func,  sizeof(double) * nw_pts);
            phi_arr = (double *)realloc(phi_arr, sizeof(double) * nw_pts);
            x_arr   = (double *)realloc(x_arr,   sizeof(double) * nw_pts);

            /*  Recompute rho, phi, and the window function.                  */
            for (j=0; j<nw_pts; ++j)
            {
                x_arr[j]   = dlp->rho_km_vals[center+j-(nw_pts-1)/2];
                phi_arr[j] = dlp->phi_rad_vals[center+j-(nw_pts-1)/2];
                w_func[j]  = fw(x_arr[j] - dlp->rho_km_vals[center], w_init);
            }
        }
        else 
        {
            /*  Adjust rho and phi to the new range.                          */
            for (j=0; j<nw_pts; ++j)
            {
                x_arr[j]   = dlp->rho_km_vals[center+j-(nw_pts-1)/2];
                phi_arr[j] = dlp->phi_rad_vals[center+j-(nw_pts-1)/2];
            }
        }

        /*  Compute the fresnel tranform about the current point.             */
        dlp->T_out[i] = FresT(x_arr, phi_arr, dlp->T_in, w_func,
                              dlp->kd_vals[center], dlp->rho_km_vals[center],
                              dlp->B_rad_vals[center], dlp->D_km_vals[center],
                              EPS, toler, dx, dlp->F_km_vals[center],
                              nw_pts, center);

        /*  Increment pointers using pointer arithmetic.                      */
        center += 1;
    }

    /*  Free variables allocated by malloc.                                   */
    free(x_arr);
    free(phi_arr);
    free(w_func);
}

void DiffractionCorrectionPerturbedNewton(DLPObj *dlp)
{
    /*  If everything executes smoothly, status should remain at zero.        */
    dlp->status = 0;

    /*  Check that the pointers to the data are not NULL.                     */
    if (!(check_dlp_data(dlp)))
    {
        /*  One of the variables has null data, return to calling function.   */
        dlp->status = 1;
        return;
    }

    /*  Variables for indexing, nw_pts is the number of points in the window. */
    long i, j, nw_pts, center;

    /*  Maximum number of iterations allowed in Newton-Raphson.               */
    long toler;

    /*  Some more necessary variables.                                        */
    double w_init, dx, two_dx;

    /*  The maximum allowed error in Newton-Raphson.                          */
    double EPS;

    /*  Set toler to 5 and EPS to e-4, both reasonable for almost all cases.  */
    toler = 5;
    EPS = 1.E-4;

    /*  Function pointers for the window function and Fresnel transform.      */
    double (*fw)(double, double);
    complex double (*FresT)(double *, double *, complex double *, double *,
                            double, double, double, double, double, long,
                            double, double, long, long, double [5]);

    /*  Cast the selected window type to the fw pointer.                      */
    if      (dlp->wtype == 0) fw = &Rect_Window_Double;
    else if (dlp->wtype == 1) fw = &Coss_Window_Double;
    else if (dlp->wtype == 2) fw = &Kaiser_Bessel_2_0_Double;
    else if (dlp->wtype == 3) fw = &Kaiser_Bessel_2_5_Double;
    else if (dlp->wtype == 4) fw = &Kaiser_Bessel_3_5_Double;
    else if (dlp->wtype == 5) fw = &Modified_Kaiser_Bessel_2_0_Double;
    else if (dlp->wtype == 6) fw = &Modified_Kaiser_Bessel_2_5_Double;
    else                      fw = &Modified_Kaiser_Bessel_3_5_Double;

    /*  Select the correct transform function.                                */
    if (dlp->use_norm) FresT = &Fresnel_Transform_Perturbed_Newton_Norm_Double;
    else               FresT = &Fresnel_Transform_Perturbed_Newton_Double;

    /* Compute first window width and window function.                        */
    center = dlp->start;

    /*  If forward calculation is selected, negate kd_vals. This has the same *
     *  effect as computation the forward model computation.                  */
    if (dlp->use_fwd)
    {
        /*  Loop over kd_vals and negate each term.                           */
        for (i=0; i<=dlp->n_used; ++i)
        {
            dlp->kd_vals[center+i] *= -1.0;
        }
    }

    /*  Compute some more necessary data.                                     */
    w_init  = dlp->w_km_vals[center];
    dx      = dlp->rho_km_vals[center+1] - dlp->rho_km_vals[center];
    two_dx  = 2.0*dx;
    nw_pts  = 2*((long)(w_init / (2.0 * dx)))+1;

    if (center - ((nw_pts-1)/2) < 0)
    {
        /*  Invalid index for the for the first point. That is, the window    *
         *  of integration required for the first point goes beyond the       *
         *  available data. Return to calling function.                       */
        dlp->n_used = (nw_pts-1)/2;
        dlp->status = 2;
        return;
    }

    /*  Allocate memory for rho, phi, and the window function.                */
    double *x_arr   = (double *)malloc(sizeof(double) * nw_pts);
    double *phi_arr = (double *)malloc(sizeof(double) * nw_pts);
    double *w_func  = (double *)malloc(sizeof(double) * nw_pts);

    /*  Check that malloc was successfull.                                    */
    if (!(x_arr) || !(w_func) || !(phi_arr))
    {
        /*  Malloc failed, return to calling function.                        */
        dlp->status = 3;
        return;
    }

    /*  Compute rho, phi, and the window function.                            */
    for (j=0; j<nw_pts; ++j)
    {
        x_arr[j]   = dlp->rho_km_vals[center+j-(nw_pts-1)/2];
        phi_arr[j] = dlp->phi_rad_vals[center+j-(nw_pts-1)/2];
        w_func[j]  = fw(x_arr[j] - dlp->rho_km_vals[center], w_init);
    }

    for (i=0; i<=dlp->n_used; ++i)
    {
        /*  If the window width changes significantly, recompute w_func.      */
        if (fabs(w_init - dlp->w_km_vals[center]) >= two_dx)
        {
            /* Reset w_init and recompute window function.                    */
            w_init  = dlp->w_km_vals[center];
            nw_pts  = 2*((int)(w_init / (2.0 * dx)))+1;

            /*  Reallocate memory since the sizes changed.                    */
            w_func  = (double *)realloc(w_func,  sizeof(double) * nw_pts);
            phi_arr = (double *)realloc(phi_arr, sizeof(double) * nw_pts);
            x_arr   = (double *)realloc(x_arr,   sizeof(double) * nw_pts);

            /*  Recompute rho, phi, and the window function.                  */
            for (j=0; j<nw_pts; ++j)
            {
                x_arr[j]   = dlp->rho_km_vals[center+j-(nw_pts-1)/2];
                phi_arr[j] = dlp->phi_rad_vals[center+j-(nw_pts-1)/2];
                w_func[j]  = fw(x_arr[j] - dlp->rho_km_vals[center], w_init);
            }
        }
        else
        {
            /*  Compute rho and phi about the new range.                      */
            for (j=0; j<nw_pts; ++j)
            {
                x_arr[j]   = dlp->rho_km_vals[center+j-(nw_pts-1)/2];
                phi_arr[j] = dlp->phi_rad_vals[center+j-(nw_pts-1)/2];
            }
        }

        /*  Compute the fresnel tranform about the current point.             */
        dlp->T_out[i] = FresT(x_arr, phi_arr, dlp->T_in, w_func,
                              dlp->kd_vals[center], dlp->rho_km_vals[center],
                              dlp->B_rad_vals[center], dlp->D_km_vals[center],
                              EPS, toler, dx, dlp->F_km_vals[center],
                              nw_pts, center, dlp->perturb);

        /*  Increment pointers using pointer arithmetic.                      */
        center += 1;
    }

    /*  Free variables allocated by malloc.                                   */
    free(x_arr);
    free(phi_arr);
    free(w_func);
}

void DiffractionCorrectionEllipse(DLPObj *dlp)
{
    /*  If everything executes smoothly, status should remain at zero.        */
    dlp->status = 0;

    /*  Check that the pointers to the data are not NULL.                     */
    if (!(check_dlp_data(dlp)))
    {
        /*  One of the variables has null data, return to calling function.   */
        dlp->status = 1;
        return;
    }

    long i, j, nw_pts, toler, center;
    double w_init, dx, two_dx, EPS;

    toler = 5;
    EPS = 1.E-4;

    double (*fw)(double, double);
    complex double (*FresT)(double *, double *, complex double *, double *,
                            double, double, double, double, double, long,
                            double, double, long, long, double, double);

    /*  Cast the selected window type to the fw pointer.                      */
    if      (dlp->wtype == 0) fw = &Rect_Window_Double;
    else if (dlp->wtype == 1) fw = &Coss_Window_Double;
    else if (dlp->wtype == 2) fw = &Kaiser_Bessel_2_0_Double;
    else if (dlp->wtype == 3) fw = &Kaiser_Bessel_2_5_Double;
    else if (dlp->wtype == 4) fw = &Kaiser_Bessel_3_5_Double;
    else if (dlp->wtype == 5) fw = &Modified_Kaiser_Bessel_2_0_Double;
    else if (dlp->wtype == 6) fw = &Modified_Kaiser_Bessel_2_5_Double;
    else                      fw = &Modified_Kaiser_Bessel_3_5_Double;

    if (dlp->use_norm) FresT = &Fresnel_Transform_Ellipse_Norm_Double;
    else               FresT = &Fresnel_Transform_Ellipse_Double;

    /* Compute first window width and window function. */
    center = dlp->start;

    if (dlp->use_fwd){
        for (i=0; i<=dlp->n_used; ++i){
            dlp->kd_vals[center+i] *= -1.0;
        }
    }

    w_init  = dlp->w_km_vals[center];
    dx      = dlp->rho_km_vals[center+1] - dlp->rho_km_vals[center];
    two_dx  = 2.0*dx;
    nw_pts  = 2*((long)(w_init / (2.0 * dx)))+1;

    if (center - ((nw_pts-1)/2) < 0){

        /*  Invalid index for the for the first point. That is, the window    *
         *  of integration required for the first point goes beyond the       *
         *  available data. Return to calling function.                       */
        dlp->n_used = (nw_pts-1)/2;
        dlp->status = 2;
        return;
    }

    double *x_arr   = (double *)malloc(sizeof(double) * nw_pts);
    double *phi_arr = (double *)malloc(sizeof(double) * nw_pts);
    double *w_func  = (double *)malloc(sizeof(double) * nw_pts);

    /*  Check that malloc was successfull.                                    */
    if (!(x_arr) || !(w_func) || !(phi_arr)){

        /*  Malloc failed, return to calling function.                        */
        dlp->status = 3;
        return;
    }

    for (j=0; j<nw_pts; ++j){
        x_arr[j]   = dlp->rho_km_vals[center+j-(nw_pts-1)/2];
        phi_arr[j] = dlp->phi_rad_vals[center+j-(nw_pts-1)/2];
        w_func[j]  = fw(x_arr[j] - dlp->rho_km_vals[center], w_init);
    }

    for (i=0; i<=dlp->n_used; ++i){

        /*  If the window width changes significantly, recompute w_func.  */
        if (fabs(w_init - dlp->w_km_vals[center]) >= two_dx) {
            // Reset w_init and recompute window function.
            w_init  = dlp->w_km_vals[center];
            nw_pts  = 2*((int)(w_init / (2.0 * dx)))+1;
            w_func  = (double *)realloc(w_func,  sizeof(double) * nw_pts);
            phi_arr = (double *)realloc(phi_arr, sizeof(double) * nw_pts);
            x_arr   = (double *)realloc(x_arr,   sizeof(double) * nw_pts);
            for (j=0; j<nw_pts; ++j){
                x_arr[j]   = dlp->rho_km_vals[center+j-(nw_pts-1)/2];
                phi_arr[j] = dlp->phi_rad_vals[center+j-(nw_pts-1)/2];
                w_func[j]  = fw(x_arr[j] - dlp->rho_km_vals[center], w_init);
            }
        }
        else {
            for (j=0; j<nw_pts; ++j){
                x_arr[j]   = dlp->rho_km_vals[center+j-(nw_pts-1)/2];
                phi_arr[j] = dlp->phi_rad_vals[center+j-(nw_pts-1)/2];
            }
        }

        /*  Compute the fresnel tranform about the current point.   */
        dlp->T_out[i] = FresT(x_arr, phi_arr, dlp->T_in, w_func,
                              dlp->kd_vals[center], dlp->rho_km_vals[center],
                              dlp->B_rad_vals[center], dlp->D_km_vals[center],
                              EPS, toler, dx, dlp->F_km_vals[center],
                              nw_pts, center, dlp->ecc, dlp->peri);

        /*  Increment pointers using pointer arithmetic.                      */
        center += 1;
    }
    free(x_arr);
    free(phi_arr);
    free(w_func);
}