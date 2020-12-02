/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

rssringoccs_ComplexDouble
rssringoccs_CDouble_Erfc(rssringoccs_ComplexDouble z)
{
    double z_x, z_y, w_x, w_y;
    double mRe_z2, mIm_z2;
    rssringoccs_ComplexDouble w, temp;

    z_x = rssringoccs_CDouble_Real_Part(z);
    z_y = rssringoccs_CDouble_Imag_Part(z);

    if (z_x == 0.0)
    {
        w_x = 1.0;

        /*  handle y -> Inf limit manually, since exp(y^2) -> Inf but         *
         *  Im[w(y)] -> 0, so IEEE will give us a NaN when it should be Inf   */
        if (z_y*z_y > 720)
        {
            if (z_y>0)
                w_y = -rssringoccs_Infinity;
            else
                w_y = rssringoccs_Infinity;
        }
        else
            w_y = -exp(z_y*z_y) * rssringoccs_Double_Faddeeva_Im(z_y);

        w = rssringoccs_CDouble_Rect(w_x, w_y);
    }
    else if (z_y == 0.0)
    {
        /*  Underflow.                                                        */
        if (z_x*z_x > 750)
        {
            w_y = -z_y;
            if (z_x >= 0.0)
                w_x = 0.0;
            else
                w_x = 2.0;

            w = rssringoccs_CDouble_Rect(w_x, w_y);
        }
        else
        {
            w_y = -z_y;
            if (z_x >= 0.0)
                w_x = exp(-z_x*z_x) * rssringoccs_Double_Erfcx(z_x);
            else
                w_x = 2.0 - exp(-z_x*z_x) * rssringoccs_Double_Erfcx(z_x);

            w = rssringoccs_CDouble_Rect(w_x, w_y);
        }
    }
    else
    {
        mRe_z2 = (z_y - z_x) * (z_x + z_y);

        /*  Im(-z^2).                                                         */
        mIm_z2 = -2*z_x*z_y;

        /*  Underflow case.                                                   */
        if (mRe_z2 < -750)
        {
            w_y = 0.0;
            if (z_x >= 0.0)
                w_x = 0.0;
            else
                w_x = 2.0;
            w = rssringoccs_CDouble_Rect(w_x, w_y);
        }
        else
        {
            if (z_x >= 0.0)
            {
                temp = rssringoccs_CDouble_Rect(mRe_z2, mIm_z2);
                temp = rssringoccs_CDouble_Exp(temp);
                w = rssringoccs_CDouble_Rect(-z_y, z_x);
                w = rssringoccs_CDouble_Faddeeva(w);
                w = rssringoccs_CDouble_Multiply(w, temp);
            }
            else
            {
                temp = rssringoccs_CDouble_Rect(mRe_z2, mIm_z2);
                temp = rssringoccs_CDouble_Exp(temp);
                w = rssringoccs_CDouble_Rect(z_y, -z_x);
                w = rssringoccs_CDouble_Faddeeva(w);
                w = rssringoccs_CDouble_Multiply(w, temp);
                w = rssringoccs_CDouble_Subtract_Real(2.0, w);
            }
        }
    }
    return w;

}
