/*  Booleans provided by this library.                                        */
#include <libtmpl/include/tmpl.h>

/*  Header file with the Tau definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  Sets the default values for a Tau objects.                                */
void rssringoccs_Tau_Set_Default_Values(rssringoccs_TAUObj* tau)
{
    if (!tau)
        return;

    if (tau->error_occurred)
        return;

    /*  The Allen Deviation for the Cassini spacecraft.                       */
    tau->sigma = 2.0E-13;

    /*  Saturn geometry is assumed perfectly cylindrical for most             *
     *  computations. This is very accurate for most resolutions.             */
    tau->ecc = 0.0;
    tau->peri = 0.0;

    /*  Default resolution for the PDS is 1 kilometer.                        */
    tau->res = 1.0;

    /*  By default no perturbation is introduced to the Fresnel kernel.       */
    tau->perturb[0] = 0.0;
    tau->perturb[1] = 0.0;
    tau->perturb[2] = 0.0;
    tau->perturb[3] = 0.0;
    tau->perturb[4] = 0.0;

    /*  Default range is "all". For Saturn this is about 70,000 to 140,000 km.*
     *  For Uranus this is different. To include all possible data sets for   *
     *  all planets we set this to 1 to 400,000. Later functions will then    *
     *  reset these values to the minimum and maximum allowed in the data set.*/
    tau->rng_req[0] = 1.0;
    tau->rng_req[1] = 4.0E5;

    /*  Default epsilon precision for the Newton-Raphson method of finding    *
     *  the stationary azimuthal angle for the Fresnel kernel. Setting this   *
     *  to a larger value may result in poor reconstructions for the most     *
     *  extreme geometries (like Rev133 for Cassini data). Setting it lower   *
     *  doesn't yield much of an improvement either.                          */
    tau->EPS = 1.0E-6;

    /*  Maximum number of iterations allowed in the Newton-Raphson method.    *
     *  The MTR paper says 4 iterations is enough for the Voyager data. For   *
     *  Cassini the geometry can be a little more extreme and so 6 iterations *
     *  may be required (like Rev133, for example). We set the max to 6.      */
    tau->toler = 6U;

    /*  The default window is the Modified Kaiser-Bessel window with alpha    *
     *  parameter set to 2 pi. The modification makes it so that the edge of  *
     *  windows is precisely zero, instead of 1 / I_0(2 pi) ~ 0.01.           */
    tau->window_func = tmpl_Double_Modified_Kaiser_Bessel_2_0;

    /*  The normalized equivalent width of the selected window function.      */
    tau->normeq = KBMD20NormEQ;

    /*  Default reconstruction method: Quartic interpolation of the           *
     *  Newton-Raphson method with dpsi/dphi perturbation taken into account. */
    tau->psinum = rssringoccs_DR_NewtonDPhiQuartic;

    /*  The order is only needed if the Legendre polynomial method is chosen. *
     *  Set this to zero.                                                     */
    tau->order = 0U;

    /*  Boolean for normalizing the reconstruction by the width of the window.*/
    tau->use_norm = tmpl_True;

    /*  Boolean for performing a forward modeling computation from the        *
     *  reconstructed data. If the reconstruction is good, and if the         *
     *  resolution is decent enough, the forward model should be very close   *
     *  to the input diffracted data. This doubles the length of the          *
     *  computation, however. Use it to verify your results. For the sake of  *
     *  time, default is to not perform the forward computation.              */
    tau->use_fwd = tmpl_False;

    /*  The resolution is given as a function of the Fresnel scale, window    *
     *  width, and Allen deviation. The Allen deviation term can be inverted  *
     *  in terms of the Lambert W function (provided by libtmpl). This        *
     *  inversion factor is called the "b-factor". Computing it makes         *
     *  the computation ever-so-slightly slower. For Cassini it is not needed,*
     *  the resulting computation is nearly the same if the b-factor is       *
     *  skipped. For Voyager it is needed. Since the time difference really   *
     *  is negligible (the reconstruction takes up about 99.9% of the         *
     *  computational time, computing window width is essentially instant)    *
     *  the default is to use the b-factor.                                   */
    tau->bfac = tmpl_True;

    /*  Several functions will print messages throughout the computation if   *
     *  the verbose Boolean is set to True. Default is silent, set to False.  */
    tau->verbose = tmpl_False;

    /*  Boolean for keeping track of errors. This starts as false. Every      *
     *  function that takes in a Tau object will check if this is True and    *
     *  abort the computation if so.                                          */
    tau->error_occurred = tmpl_False;

    /*  An error message is stored in the event of an error to let the user   *
     *  know where the computation went wrong. Default is a NULL pointer for  *
     *  no message at all.                                                    */
    tau->error_message = NULL;

    /*  These values are determined by the DLP data. Set them to their zero   *
     *  values just to have them initialized.                                 */
    tau->dx_km = 0.0;
    tau->rng_list[0] = 0.0;
    tau->rng_list[1] = 0.0;
}
/*  End of rssringoccs_Tau_Set_Default_Values.                                */
