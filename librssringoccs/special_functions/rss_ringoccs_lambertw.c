/*  The C Standard Library header for math functions and more found here.     */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

float rssringoccs_Float_LambertW(float x)
{
    float EPS = 1.0e-4F;
    float x0, dx, exp_x0;
    if ((x < rssringoccs_Infinity_F) && (x>-rssringoccs_Rcpr_Euler_E))
    {
        if (x > 2.0F)
            x0 = rssringoccs_Float_Log(x/rssringoccs_Float_Log(x));
        else if (x > -rssringoccs_Rcpr_Euler_E)
            x0 = x;
        else if (x == -rssringoccs_Rcpr_Euler_E)
            return -1.0F;
        else
            return rssringoccs_NaN_F;

        exp_x0 = rssringoccs_Float_Exp(x0);
        dx = (x0*exp_x0-x)/(exp_x0*(x0+1.0) -
                            (x0+2.0)*(x0*exp_x0-x)/(2.0*x0+2.0));
        x0 = x0 - dx;

        while (rssringoccs_Float_Abs(dx) > EPS)
        {
            exp_x0 = rssringoccs_Float_Exp(x0);
            dx = (x0*exp_x0-x)/(exp_x0*(x0+1.0) -
                                (x0+2.0)*(x0*exp_x0-x)/(2.0*x0+2.0));
            x0 = x0 - dx;
        }
        return x0;
    }
    else if (x==-rssringoccs_Rcpr_Euler_E)
        return -1.0F;
    else if (x<-rssringoccs_Rcpr_Euler_E)
        return rssringoccs_NaN;
    else
        return rssringoccs_Infinity;
}

double rssringoccs_Double_LambertW(double x)
{
    double EPS = 1.0e-8;
    double x0, dx, exp_x0;
    if ((x < rssringoccs_Infinity) && (x>-rssringoccs_Rcpr_Euler_E))
    {
        if (x > 2.0){
            x0 = rssringoccs_Double_Log(x/rssringoccs_Double_Log(x));
        }
        else if (x > -rssringoccs_Rcpr_Euler_E){
            x0 = x;
        }
        else if (x == -rssringoccs_Rcpr_Euler_E){
            return -1.0;
        }
        else {
            return rssringoccs_NaN;
        }
        exp_x0 = rssringoccs_Double_Exp(x0);
        dx = (x0*exp_x0-x)/(exp_x0*(x0+1.0) -
                            (x0+2.0)*(x0*exp_x0-x)/(2.0*x0+2.0));
        x0 = x0 - dx;
        while (rssringoccs_Double_Abs(dx) > EPS)
        {
            exp_x0 = rssringoccs_Double_Exp(x0);
            dx = (x0*exp_x0-x)/(exp_x0*(x0+1.0) -
                                (x0+2.0)*(x0*exp_x0-x)/(2.0*x0+2.0));
            x0 = x0 - dx;
        }
        return x0;
    }
    else if (x==-rssringoccs_Rcpr_Euler_E)
        return -1.0;
    else if (x<-rssringoccs_Rcpr_Euler_E)
        return rssringoccs_NaN;
    else
        return rssringoccs_Infinity;
}

long double rssringoccs_LDouble_LambertW(long double x)
{
    long double EPS = 1.0e-16;
    long double x0, dx, exp_x0;
    if ((x < rssringoccs_Infinity_L) && (x>-rssringoccs_Rcpr_Euler_E))
    {
        if (x > 2.0){
            x0 = rssringoccs_LDouble_Log(x/rssringoccs_LDouble_Log(x));
        }
        else if (x > -rssringoccs_Rcpr_Euler_E){
            x0 = x;
        }
        else if (x == -rssringoccs_Rcpr_Euler_E){
            return -1.0;
        }
        else {
            return rssringoccs_NaN_L;
        }
        exp_x0 = rssringoccs_LDouble_Exp(x0);
        dx = (x0*exp_x0-x)/(exp_x0*(x0+1.0) -
                            (x0+2.0)*(x0*exp_x0-x)/(2.0*x0+2.0));
        x0 = x0 - dx;
        while (rssringoccs_LDouble_Abs(dx) > EPS){
            exp_x0 = rssringoccs_LDouble_Exp(x0);
            dx = (x0*exp_x0-x)/(exp_x0*(x0+1.0) -
                                (x0+2.0)*(x0*exp_x0-x)/(2.0*x0+2.0));
            x0 = x0 - dx;
        }
        return x0;
    }
    else if (x==-rssringoccs_Rcpr_Euler_E)
        return -1.0;
    else if (x<-rssringoccs_Rcpr_Euler_E)
        return rssringoccs_NaN_L;
    else
        return rssringoccs_Infinity_L;
}
