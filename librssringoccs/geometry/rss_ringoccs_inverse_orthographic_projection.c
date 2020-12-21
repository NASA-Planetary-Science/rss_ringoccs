#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

rssringoccs_ThreeVector
rssringoccs_Inverse_Orthographic_Projection(rssringoccs_TwoVector P,
                                            rssringoccs_ThreeVector u)
{
    /*  Declare all necessary variables. C89 requires this at the top.        */
    double x, y, z, radius;
    rssringoccs_ThreeVector out, X, Y, u_hat;

    /*  Extract the X and Y components from the point P.                      */
    x = rssringoccs_TwoVector_X(P);
    y = rssringoccs_TwoVector_Y(P);

    /*  The radius of the sphere we'll be computing with is just the norm     *
     *  of the input ThreeVector u, so compute this.                          */
    radius = rssringoccs_ThreeVector_Euclidean_Norm(u);

    /*  If the norm of P is greater than the radius the inverse stereographic *
     *  projection is undefined. We'll return Not-a-Number in this case.      */
    if ((x*x + y*y) > radius)
    {
        out = rssringoccs_ThreeVector_Rect(rssringoccs_NaN,
                                           rssringoccs_NaN,
                                           rssringoccs_NaN);
    }
    else
    {
        /*  Normalize the input u vector so that it lies on the sphere.       */
        u_hat = rssringoccs_ThreeVector_Normalize(u);

        /*  Get a vector orthogonal to u and normalize it.                    */
        X = rssringoccs_Orthogonal_ThreeVector(u);
        X = rssringoccs_ThreeVector_Normalize(X);

        /*  Compute the cross product of X and u, giving as an orthonormal    *
         *  basis of three dimensional space: (X, Y, u_hat).                  */
        Y = rssringoccs_Cross_Product(X, u_hat);

        /*  The z component of our sphere is chosen so that x^2+y^2+z^2=r^2   *
         *  and so that it is positive.                                       */
        z = rssringoccs_Double_Sqrt(radius*radius - x*x - y*y);

        /*  The point on the sphere now satisfies x*X + y*Y + z*u_hat. We     *
         *  compute this and return.                                          */
        out = rssringoccs_ThreeVector_Add(
            rssringoccs_ThreeVector_Add(
                rssringoccs_ThreeVector_Scale(x, X),
                rssringoccs_ThreeVector_Scale(y, Y)
            ),
            rssringoccs_ThreeVector_Scale(z, u_hat)
        );
    }

    return out;
}
