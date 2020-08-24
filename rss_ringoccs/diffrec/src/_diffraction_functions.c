/******************************************************************************
 *                                 LICENSE                                    *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify it   *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  rss_ringoccs is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************
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
 *          Simple FFT Method:                                                *
 *              Uses an FFT to compute the Fresnel transform. This assumes    *
 *              the geometry about the center of the data set is an accurate  *
 *              representative of the entire occultation. Hence, this method  *
 *              accurate around the center and inaccurate towards the edges.  *
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
 *                                                                            *
 *      For near-constant geometry, or regions where the geometry does not    *
 *      depend to heavily on phi, psi is approximately a function of the form *
 *      psi(rho,rho_0) = f(rho-rho_0) and so this integral is a convolution   *
 *      of w * exp(-i psi) and T_hat. Using the convolution theorem, we get:  *
 *                                                                            *
 *      F(T)   =    F(T_hat) * F(w * exp(-i psi))                             *
 *                                                                            *
 *      T      =    F^-1(F(T_hat) * F(w * exp(-i psi)))                       *
 *                                                                            *
 *      Where F(f) is the Fourier transform of f, and F^-1 is the inverse     *
 *      Fourier transform.                                                    *
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
 *          Norm =  _________________                        *
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
 *                                                                            *
 *  In addition, for full use of all the routines one will need the FFTW      *
 *  library. Installation instruction can be found in the rss_ringoccs PDF.   *
 ******************************************************************************
 *  It is anticipated that many users of this code will have experience in    *
 *  either Python or IDL, but not C. Many comments are left to explain as     *
 *  much as possible. Vagueness or unclear code should be reported to:        *
 *  https://github.com/NASA-Planetary-Science/rss_ringoccs/issues             *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       June 21, 2019                                                 *
 ******************************************************************************
 *                                History                                     *
 ******************************************************************************
 *  2019/06/21 (Ryan Maguire):                                                *
 *      Initial commit.                                                       *
 *  2020/07/28 (Ryan Maguire):                                                *
 *      Clarified comments, fixed bug in error checks.                        *
 *  2020/08/22 (Ryan Maguire):                                                *
 *      Added FFT routine.                                                    *
 ******************************************************************************/

/*  The malloc and realloc functions are contained in stdlib.h. Various       *
 *  diffraction-based functions are found in _fresnel_kernel.h and            *
 *  _window_functions.h. The DLPObj is defined in _diffraction_functions.h.   *
 *  This file includes <complex.h> which must be include BEFORE <fftw3.h>.    *
 *  Do NOT swap the order of these include files.                             */
#include "_diffraction_functions.h"

/*  All of the kernel functions (psi) found here.                             */
#include "_fresnel_kernel.h"

/*  C Library for using FFT routines. This is NOT part of the standard C      *
 *  and you'll need to compile/build before compiling/using these routines.   */
#include <fftw3.h>

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
 *      fw (double (*)(double, double)):                                      *
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

/******************************************************************************
 *  Function:                                                                 *
 *      check_dlp_data                                                        *
 *  Purpose:                                                                  *
 *      Check if the input DLP data is valid.                                 *
 *  Arguments:                                                                *
 *      dlp (DLPObj *):                                                       *
 *          An instance of the DLPObj structure defined in                    *
 *          _diffraction_functions.h.                                        *
 *  Output:                                                                   *
 *      Out (Boolean):                                                        *
 *          An integer, either 0 or 1 depending on if the checks passed.      *
 *  Notes:                                                                    *
 *      1.) This is a "static" function meanings its use is limited to this   *
 *          file. Its only real purpose is for checking DLP instances which   *
 *          are the inputs to the various diffraction correction functions.   *
 ******************************************************************************/
static int check_dlp_data(DLPObj *dlp)
{
    /*  If any of these pointers are void, return 0. Else, pass.              */
    if ((!dlp->T_in) || (!dlp->rho_km_vals)  || (!dlp->F_km_vals)
                     || (!dlp->phi_rad_vals) || (!dlp->kd_vals)
                     || (!dlp->B_rad_vals)   || (!dlp->D_km_vals)
                     || (!dlp->w_km_vals)    || (!dlp->T_out)) return 0;
    else return 1;
}

/******************************************************************************
 *  Function:                                                                 *
 *      check_data_range                                                      *
 *  Purpose:                                                                  *
 *      Check if the window ranges needed for reconstruction are permitted by *
 *      the available data. This check is important to avoid segmentation     *
 *      faults at runtime.                                                    *
 *  Arguments:                                                                *
 *      dlp (DLPObj *):                                                       *
 *          An instance of the DLPObj structure defined in                    *
 *          _diffraction_functions.h.                                        *
 *  Output:                                                                   *
 *      Out (Boolean):                                                        *
 *          An integer, either 0 or 1 depending on if the checks passed.      *
 *  Notes:                                                                    *
 *      1.) This is a "static" function meanings its use is limited to this   *
 *          file. Its only real purpose is for checking DLP instances which   *
 *          are the inputs to the various diffraction correction functions.   *
 *      2.) If the check fails, there is something wrong with the input data. *
 *          Removing this check may result in a "segmentation fault 11" error *
 *          at runtime with no additional error message, a scenario one would *
 *          like to avoid.                                                    *
 ******************************************************************************/
static int check_data_range(DLPObj *dlp, double two_dx)
{
    /* Create variables to check what the minimum and maximum indices are.    */
    long current_min, current_max;

    /* Variable for keeping track of the largest window needed in the loop.   */
    long max_nw_pts = 0;

    /* The index corresponding to the largest window.                         */
    long worst_point_min = dlp->start;
    long worst_point_max = dlp->start;

    /* Variables to keep track of the min and max window indices.             */
    long min_requested = dlp->start;
    long max_requested = min_requested + dlp->n_used;

    /* The number of elements to check.                                       */
    long N = max_requested;

    /* Variable for indexing the for loop.                                    */
    long i;

    /* Variables for the window size and number of points in the window.      */
    double win_size;
    long nw_pts;

    /* Loop through every point, check window width, and ensure you have      *
     * enough data to the left and right for data processing.                 */
    for (i=dlp->start; i<=N; ++i){
        /*  Gather the window size and number of points in the window.        */
        win_size = dlp->w_km_vals[i];
        nw_pts = ((long)(win_size / two_dx))+1;

        /*  The smallest radius needed for the window corresponds to the      *
         *  current point minus the number of points in the window. Similarly *
         *  the maximum corresponds to the sum. Compute these.                */
        current_min = i-nw_pts;
        current_max = i+nw_pts;

        /*  If current_min is smaller than min_requested reset min_requrested *
         *  and similarly if current_max is larger than max_requested.        */
        if (current_min < min_requested){
            min_requested = current_min;
            worst_point_min = i;
        }
        else if (current_max < max_requested){
            worst_point_max = i;
            max_requested = current_max;
        }
        if (nw_pts > max_nw_pts) max_nw_pts = nw_pts;
    }

    /*  If min_requested is negative, the window is too larger. Similarly, if *
     *  max_requested goes between the size of the array.                     */
    if (min_requested < 0){
        dlp->start  = worst_point_min;
        dlp->n_used = max_nw_pts;
        return 0;
    }
    else if (max_requested > dlp->arr_size){
        dlp->start  = worst_point_max;
        dlp->n_used = max_nw_pts;
        return 0;
    }
    else return 1;
}

/******************************************************************************
 *  Function:                                                                 *
 *      DiffractionCorrectionFresnel                                          *
 *  Purpose:                                                                  *
 *      Compute the Fresnel transform using the classic Fresnel quadratic     *
 *      approximation to the fresnel kernel.                                  *
 *  Arguments:                                                                *
 *      dlp (DLPObj *):                                                       *
 *          An instance of the DLPObj structure defined in                    *
 *          _diffraction_correction.h. This contains all of the necessary    *
 *          data for diffraction correction, including the geometry of the    *
 *          occultation and actual power and phase data.                      *
 *  Output:                                                                   *
 *      Nothing (void):                                                       *
 *          This is a void function, so no actual output is provided. However *
 *          the T_out pointer within the dlp structure will be changed at the *
 *          end, containing the diffraction correction data.                  *
 *  Notes:                                                                    *
 *      1.) This code uses the Fresnel approximation which has been known to  *
 *          fail for several different occultations, especially ones of very  *
 *          low angle (small B values). Take this into consideration when     *
 *          performing any analysis.                                          *
 *      2.) While this may be inaccurate for certain occultations, it is      *
 *          immensely fast, capable of processing the entire Rev007 E         *
 *          in less than a second at 1km resolution.                          *
 ******************************************************************************/
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

    /*  Retrieve the starting point from the DLPObj instance.                 */
    center = dlp->start;

    /*  Compute some extra necessary variables.                               */
    w_init  = dlp->w_km_vals[center];
    dx      = dlp->rho_km_vals[center+1] - dlp->rho_km_vals[center];
    two_dx  = 2.0*dx;
    nw_pts  = ((long)(w_init / two_dx))+1;

    // Check to ensure you have enough data to the left.
    if (!check_data_range(dlp, two_dx))
    {
        /*  One of the points has too large of a window width to process.     *
         *  Returning with error message.                                     */
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

        /*  If forward transform is selected, negate x_arr. This is a quick   *
         *  way to multiply by +/- 1. (dlp->use_fwd == 0) is a Boolean which  *
         *  will be 1 if the statement is true and 0 if false. Similarly for  *
         *  (dlp->use_fwd == 1). Hence this returns +/- 1.                    */
        x_arr[j] *= (dlp->use_fwd == 0) - (dlp->use_fwd == 1);
    }

    /*  Compute the Fresnel transform across the input data.                  */
    for (i=0; i<=dlp->n_used; ++i)
    {
        /*  If the window width has deviated more the 2*dx, reset variables.  *
         *  fabs is the absolute value function for double precision          *
         *  variables and is defined in the built-in math.h.                  */
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

/******************************************************************************
 *  Function:                                                                 *
 *      DiffractionCorrectionLegendre                                         *
 *  Purpose:                                                                  *
 *      Compute the Fresnel transform using Legendre polynomials to           *
 *      approximate the fresnel kernel.                                       *
 *  Arguments:                                                                *
 *      dlp (DLPObj *):                                                       *
 *          An instance of the DLPObj structure defined in                    *
 *          _diffraction_correction.h. This contains all of the necessary    *
 *          data for diffraction correction, including the geometry of the    *
 *          occultation and actual power and phase data.                      *
 *  Output:                                                                   *
 *      Nothing (void):                                                       *
 *          This is a void function, so no actual output is provided. However *
 *          the T_out pointer within the dlp structure will be changed at the *
 *          end, containing the diffraction correction data.                  *
 *  Notes:                                                                    *
 *      1.) This routine allows for any selection of polynomial of degree     *
 *          greater than or equal to 2, though for degree 2 it is better to   *
 *          use the Fresnel option since the symmetry nearly doubles the      *
 *          speed of the computation. For anything higher than degree 8 there *
 *          is no real change in the accuracy, even for low inclination       *
 *          occultation observations.                                         *
 *      2.) Like the Fresnel approximation, the Legendre approximation has    *
 *          issues reconstructing data at low B angles. This is because the   *
 *          Legendre approximation assumes the first iteration of the Newton  *
 *          Raphson method is good enough, whereas in reality 3-4 iterations  *
 *          may be needed, like in Rev133.                                    *
 ******************************************************************************/
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
        poly_order = dlp->order;

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

    /*  If forward tranform is set, negate the kd_vals variable. This has     *
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
    w_init = dlp->w_km_vals[center];
    dx     = dlp->rho_km_vals[center+1] - dlp->rho_km_vals[center];
    two_dx = 2.0*dx;
    nw_pts = (long)(w_init / two_dx)+1;

    /* Check to ensure you have enough data to the left.                      */
    if (!check_data_range(dlp, two_dx))
    {
        /*  One of the points has too large of a window width to process.     *
         *  Returning with error message.                                     */
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
     *  that it's values range from -W/2 to zero, W being the window width.   */
    if (!(x_arr)    ||    !(w_func)            ||    !(legendre_p)
                    ||    !(alt_legendre_p)    ||    !(fresnel_ker_coeffs))
    {
        /*  Malloc failed, return to calling function.                        */
        dlp->status = 3;
        return;
    }
    else reset_window(x_arr, w_func, dx, w_init, nw_pts, fw);

    /* Loop through each point and begin the reconstruction.                  */
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

/******************************************************************************
 *  Function:                                                                 *
 *      DiffractionCorrectionNewton                                           *
 *  Purpose:                                                                  *
 *      Compute the Fresnel transform using the Newton-Raphson method to      *
 *      compute the stationary value of the Fresnel-Kernel.                   *
 *  Arguments:                                                                *
 *      dlp (DLPObj *):                                                       *
 *          An instance of the DLPObj structure defined in                    *
 *          _diffraction_correction.h. This contains all of the necessary    *
 *          data for diffraction correction, including the geometry of the    *
 *          occultation and actual power and phase data.                      *
 *  Output:                                                                   *
 *      Nothing (void):                                                       *
 *          This is a void function, so no actual output is provided. However *
 *          the T_out pointer within the dlp structure will be changed at the *
 *          end, containing the diffraction correction data.                  *
 *  Notes:                                                                    *
 *      1.) This method is the most accurate, but much slower than the        *
 *          Fresnel and Legendre options. It is accurate for every Rev of the *
 *          Cassini mission with the exception of the Rev133 occultation      *
 *          of which only the Ka band produces accurate results. For X and S  *
 *          bands one needs to use the Perturbed Newton method.               *
 *      2.) The polynomials from the MTR86 are available via the dlp.interp   *
 *          variable. These polynomials are slower and less accurate than the *
 *          normal Newton method since to compute the polynomials the         *
 *          Newton-Raphson method must be performed, and hence the            *
 *          polynomials increase the number of computations needed. The real  *
 *          use of them arises if one uses FFT methods. This routine does NOT *
 *          use FFTs, but rather ordinary integration.                        *
 ******************************************************************************/
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

    /*  If forward tranform is set, negate the kd_vals variable. This has     *
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

    /* Check to ensure you have enough data to the left.                      */
    if (!check_data_range(dlp, two_dx))
    {
        /*  One of the points has too large of a window width to process.     *
         *  Returning with error message.                                     */
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

/******************************************************************************
 *  Function:                                                                 *
 *      DiffractionCorrectionPerturbedNewton                                  *
 *  Purpose:                                                                  *
 *      Compute the Fresnel transform using the Newton-Raphson method to      *
 *      compute the stationary value of the Fresnel-Kernel. The user may then *
 *      perturb the Fresnel kernel with an arbitrary quartic polynomial.      *
 *  Arguments:                                                                *
 *      dlp (DLPObj *):                                                       *
 *          An instance of the DLPObj structure defined in                    *
 *          _diffraction_correction.h. This contains all of the necessary    *
 *          data for diffraction correction, including the geometry of the    *
 *          occultation and actual power and phase data.                      *
 *  Output:                                                                   *
 *      Nothing (void):                                                       *
 *          This is a void function, so no actual output is provided. However *
 *          the T_out pointer within the dlp structure will be changed at the *
 *          end, containing the diffraction correction data.                  *
 *  Notes:                                                                    *
 *      1.) This method can be very accurate, but requires guess and check    *
 *          work from the user since there is no magical input of polynomials *
 *          that works uniformly. If using for Saturn based data, the Encke   *
 *          gap is your friend. Start with the polynomials set to zero, which *
 *          is just the standard Newton-Raphson method. Then look to the      *
 *          Encke gap which lies between ~132900-134200 km. A good            *
 *          reconstruction should have produced a near perfect square well    *
 *          (see below). If it is not a square well, perturb the cubic and    *
 *          quartic terms until it is. The quadratic term is usually fine and *
 *          the linear term is nearly zero so you most likely do not need to  *
 *          waste time with these.                                            *
 *                                                                            *
 *          Power                                                             *
 *         |                     _________                            *
 *         |                     |               |                            *
 *         |                     |               |                            *
 *         |                     |               |                            *
 *         |                     |               |                            *
 *         |                     |               |                            *
 *         |                     |               |                            *
 *         | __________|               |__________        *
 *         |                         Encke Gap                                *
 *         |________________________________   *
 *                                                                    Radius  *
 ******************************************************************************/
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
     *  effect as computing the forward model.                                */
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

    /* Check to ensure you have enough data to the left.                      */
    if (!check_data_range(dlp, two_dx))
    {
        /*  One of the points has too large of a window width to process.     *
         *  Returning with error message.                                     */
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

/******************************************************************************
 *  Function:                                                                 *
 *      DiffractionCorrectionEllipse                                          *
 *  Purpose:                                                                  *
 *      Compute the Fresnel transform using the Newton-Raphson method to      *
 *      compute the stationary value of the Fresnel-Kernel and introducing a  *
 *      perturbation to account for eccentric rings.                          *
 *  Arguments:                                                                *
 *      dlp (DLPObj *):                                                       *
 *          An instance of the DLPObj structure defined in                    *
 *          _diffraction_correction.h. This contains all of the necessary    *
 *          data for diffraction correction, including the geometry of the    *
 *          occultation and actual power and phase data.                      *
 *  Output:                                                                   *
 *      Nothing (void):                                                       *
 *          This is a void function, so no actual output is provided. However *
 *          the T_out pointer within the dlp structure will be changed at the *
 *          end, containing the diffraction correction data.                  *
 *  Notes:                                                                    *
 *      1.) This function is only needed for Uranus data where the            *
 *          eccentricity is significant enough to effect the circular         *
 *          approximation employed by the main routines.                      *
 *      2.) As of now (July 2020) this routine is unable to perfectly         *
 *          reproduce the results from Gresh's paper. Correcting this is on   *
 *          the rss_ringoccs TODO list.                                       *
 ******************************************************************************/
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

    /*  i and j are used for indexing, nw_pts is number of points in window.  */
    long i, j, nw_pts, center;

    /*  More necessary variables.                                             */
    double w_init, dx, two_dx;

    /*  toler is the maximum number of iterations allowed in Newton-Raphson.  */
    long toler;

    /*  The maximum error allowed in the Newton-Raphson scheme.               */
    double EPS;

    /*  Set toler to 5 and EPS to e-4, both reasonable for most purposes.     */
    toler = 5;
    EPS = 1.E-4;

    /*  Function pointers for the window function and Fresnel transform.      */
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

    /*  Select the appropriate transform function.                            */
    if (dlp->use_norm) FresT = &Fresnel_Transform_Ellipse_Norm_Double;
    else               FresT = &Fresnel_Transform_Ellipse_Double;

    /* Compute first window width and window function. */
    center = dlp->start;

    /*  If use_fwd, negate kd_vals. This is equivalent to forward computation.*/
    if (dlp->use_fwd)
    {
        /*  Negate each term of kd_vals.                                      */
        for (i=0; i<=dlp->n_used; ++i)
        {
            dlp->kd_vals[center+i] *= -1.0;
        }
    }

    w_init  = dlp->w_km_vals[center];
    dx      = dlp->rho_km_vals[center+1] - dlp->rho_km_vals[center];
    two_dx  = 2.0*dx;
    nw_pts  = 2*((long)(w_init / (2.0 * dx)))+1;

    /* Check to ensure you have enough data to the left.                      */
    if (!check_data_range(dlp, two_dx))
    {
        /*  One of the points has too large of a window width to process.     *
         *  Returning with error message.                                     */
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

    /*  Set rho, phi, and the window function to their starting values.       */
    for (j=0; j<nw_pts; ++j)
    {
        x_arr[j]   = dlp->rho_km_vals[center+j-(nw_pts-1)/2];
        phi_arr[j] = dlp->phi_rad_vals[center+j-(nw_pts-1)/2];
        w_func[j]  = fw(x_arr[j] - dlp->rho_km_vals[center], w_init);
    }

    for (i=0; i<=dlp->n_used; ++i){

        /*  If the window width changes significantly, recompute w_func.      */
        if (fabs(w_init - dlp->w_km_vals[center]) >= two_dx)
        {
            /* Reset w_init and recompute window function.                    */
            w_init  = dlp->w_km_vals[center];
            nw_pts  = 2*((int)(w_init / (2.0 * dx)))+1;

            /*  Reallocate memory since sizes have changed.                   */
            w_func  = (double *)realloc(w_func,  sizeof(double) * nw_pts);
            phi_arr = (double *)realloc(phi_arr, sizeof(double) * nw_pts);
            x_arr   = (double *)realloc(x_arr,   sizeof(double) * nw_pts);

            /*  Recompute rho, phi, and window function.                      */
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
                              nw_pts, center, dlp->ecc, dlp->peri);

        /*  Increment pointers using pointer arithmetic.                      */
        center += 1;
    }

    /*  Free allocated memory.                                                */
    free(x_arr);
    free(phi_arr);
    free(w_func);
}

/******************************************************************************
 *  Function:                                                                 *
 *      DiffractionCorrectionSimpleFFT                                        *
 *  Purpose:                                                                  *
 *      Compute the Fresnel transform using an FFT about the center of the    *
 *      data. This is the fastest method, but assumes the geometry about the  *
 *      midpoint is an accurate representation of the entire occultation.
 *  Arguments:                                                                *
 *      dlp (DLPObj *):                                                       *
 *          An instance of the DLPObj structure defined in                    *
 *          _diffraction_correction.h. This contains all of the necessary    *
 *          data for diffraction correction, including the geometry of the    *
 *          occultation and actual power and phase data.                      *
 *  Output:                                                                   *
 *      Nothing (void):                                                       *
 *          This is a void function, so no actual output is provided. However *
 *          the T_out pointer within the dlp structure will be changed at the *
 *          end, containing the diffraction correction data.                  *
 *  Notes:                                                                    *
 *      1.) This method is the fast, but least accurate. It is very accurate  *
 *          near the midpoint, but assumes the geometry of this point is a    *
 *          fair representative of all of the geometry. The further one gets  *
 *          from the midpoint, the less accurate this is.                     *
 *      2.) This function uses FFTW, which is a NON-STANDARD C Library. Since *
 *          there is no FFT routine in the standard C library, FFTW has       *
 *          somewhat become the de facto standard.                            *
 ******************************************************************************/
void DiffractionCorrectionSimpleFFT(DLPObj *dlp)
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

    /*  Variable for indexing.                                                */
    long i, i_shift;

    /*  Some variables needed for reconstruction.                             */
    double dx, norm;
    complex double *ker;
    complex double *fft_ker;
    complex double *fft_in;
    complex double *fft_out;
    complex double *T_in;
    complex double *T_out;

    /*  Variable for the center of the data set.                              */
    long center = dlp->start + (dlp->n_used/2);

    /*  Window width of the midpoint.                                         */
    double window_width = dlp->w_km_vals[center];

    /*  Compute the distance between samples.                                 */
    dx = dlp->rho_km_vals[center+1] - dlp->rho_km_vals[center];

    /*  Number of points in half the window.                                  */
    long nw_pts = (long)(window_width/(2.0*dx));

    /*  Number of points in the data set, the start/end point and array size. */
    long data_size;

    data_size = dlp->n_used+2*nw_pts+1;

    /* Variables for shifting and keeping track of indexing.                  */
    long shift = data_size/2;
    long current_point;

    /*  Toler is the number of iterations allowed in Newton-Raphson.          */
    long toler;

    /*  EPS is the maximum allowed error in the Newton-Raphson scheme.        */
    double EPS;

    /*  Set toler to 5 and EPS to e-4, reasonable for almost all cases.       */
    toler = 5;
    EPS = 1.E-4;

    /*  Function pointers for the window function                             */
    double (*fw)(double, double);

    /*  Cast the selected window type to the fw pointer.                      */
    if      (dlp->wtype == 0) fw = &Rect_Window_Double;
    else if (dlp->wtype == 1) fw = &Coss_Window_Double;
    else if (dlp->wtype == 2) fw = &Kaiser_Bessel_2_0_Double;
    else if (dlp->wtype == 3) fw = &Kaiser_Bessel_2_5_Double;
    else if (dlp->wtype == 4) fw = &Kaiser_Bessel_3_5_Double;
    else if (dlp->wtype == 5) fw = &Modified_Kaiser_Bessel_2_0_Double;
    else if (dlp->wtype == 6) fw = &Modified_Kaiser_Bessel_2_5_Double;
    else                      fw = &Modified_Kaiser_Bessel_3_5_Double;

    /*  Scale factor for the FFT.                                             */
    complex double scale_factor = 0.5*dx*(1.0+_Complex_I)/((data_size+1)^2);

    /*  Allocate memory for the Fresnel kernel and other variables.           */
    ker     = (complex double *)malloc(sizeof(complex double)*(data_size));
    fft_ker = (complex double *)malloc(sizeof(complex double)*(data_size));
    fft_in  = (complex double *)malloc(sizeof(complex double)*(data_size));
    fft_out = (complex double *)malloc(sizeof(complex double)*(data_size));
    T_in    = (complex double *)malloc(sizeof(complex double)*(data_size));
    T_out   = (complex double *)malloc(sizeof(complex double)*(data_size));

    /*  The fresnel kernel.                                                   */
    double psi;

    /*  Stationary value for the Fresnel kernel.                              */
    double phi;

    /*  Independent variable for the window function.                         */
    double window_func_x;

    /* Radius of the midpoint.                                                */
    double rho = dlp->rho_km_vals[center];

    /*  Complex exponential of Fresnel kernel, weighted by window function.   */
    complex double arg;

    /*  Compute the windowing function and Psi.                               */
    for (i=0; i < data_size; ++i)
    {
        current_point = dlp->start+i-nw_pts;
        window_func_x = rho-dlp->rho_km_vals[current_point];
        phi = Newton_Raphson_Fresnel_Psi(dlp->kd_vals[current_point], rho,
                                         dlp->rho_km_vals[current_point],
                                         dlp->phi_rad_vals[current_point],
                                         dlp->phi_rad_vals[current_point],
                                         dlp->B_rad_vals[current_point],
                                         dlp->D_km_vals[current_point],
                                         EPS, toler);
        psi = -Fresnel_Psi_Double(dlp->kd_vals[current_point], rho,
                                  dlp->rho_km_vals[current_point], phi,
                                  dlp->phi_rad_vals[current_point],
                                  dlp->B_rad_vals[current_point],
                                  dlp->D_km_vals[current_point]);

        /*  If forward tranform is set, negate the Fresnel kernel.            */
        if (dlp->use_fwd) psi *= -1.0;

        arg  = cos(psi)+_Complex_I*sin(psi);
        arg *= fw(window_func_x, window_width);
        T_in[i] = dlp->T_in[current_point];
        ker[i] = arg;
    }

    /*  Select the correct Fresnel transformation.                            */
    if (dlp->use_norm)norm = 1.0;
    else norm = 1.0;

    /* Check to ensure you have enough data to the left.                      */
    if (!check_data_range(dlp, 2.0*dx))
    {
        /*  One of the points has too large of a window width to process.     *
         *  Returning with error message.                                     */
        dlp->status = 2;
        return;
    }

    /*  Check that malloc was successful.                                     */
    if (!(ker))
    {
        /*  Malloc failed, return to calling function.                        */
        dlp->status = 3;
        return;
    }

    fftw_plan p;
    p = fftw_plan_dft_1d(data_size, ker, fft_ker, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    p = fftw_plan_dft_1d(data_size, T_in, fft_in, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    for(i=0; i<=data_size; ++i) fft_out[i] = fft_ker[i]*fft_in[i];

    p = fftw_plan_dft_1d(data_size, fft_out, T_out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    for(i=0; i<=dlp->n_used; ++i)
    {
        i_shift = (nw_pts+i + shift) % (data_size);
        dlp->T_out[i] = T_out[i_shift];
        dlp->T_out[i] *= scale_factor/dlp->F_km_vals[dlp->start+i+1];
    }

    /*  Free variables allocated by malloc.                                   */
    free(ker);
    free(fft_ker);
    free(fft_in);
    free(fft_out);
    free(T_in);
    free(T_out);
}
