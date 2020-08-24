#include "special_functions.h"

complex double Square_Wave_Diffraction_Double(double x, double W, double F,
                                              long N)
{
    long i=0;

    complex double T_hat;
    double a, b;

    a = (long)(x/(2*W))*2*W;
    a -= 2*W*N;
    if (a<0) a=0;
    b = a+W;

    T_hat = Gap_Diffraction_Double(x, a, b, F);

    N *= 2;

    for (i=0; i<N; ++i){
        a += 2*W;
        b += 2*W;
        T_hat += Gap_Diffraction_Double(x, a, b, F);
    }
    return T_hat;
}
