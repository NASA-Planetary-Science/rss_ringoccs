#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

double rssringoccs_TwoVector_X(rssringoccs_TwoVector P)
{
    double x;

    x = P.dat[0];
    return x;
}

double rssringoccs_TwoVector_Y(rssringoccs_TwoVector P)
{
    double y;

    y = P.dat[1];
    return y;
}

double rssringoccs_ThreeVector_X(rssringoccs_ThreeVector P)
{
    double x;

    x = P.dat[0];
    return x;
}

double rssringoccs_ThreeVector_Y(rssringoccs_ThreeVector P)
{
    double y;

    y = P.dat[1];
    return y;
}

double rssringoccs_ThreeVector_Z(rssringoccs_ThreeVector P)
{
    double z;

    z = P.dat[2];
    return z;
}

rssringoccs_TwoVector rssringoccs_TwoVector_Rect(double x, double y)
{
    rssringoccs_TwoVector P;
    P.dat[0] = x;
    P.dat[1] = y;

    return P;
}

rssringoccs_ThreeVector
rssringoccs_ThreeVector_Rect(double x, double y, double z)
{
    rssringoccs_ThreeVector P;

    P.dat[0] = x;
    P.dat[1] = y;
    P.dat[2] = z;

    return P;
}

double rssringoccs_Euclidean_Norm_2D(rssringoccs_TwoVector P)
{
    double x, y, norm;

    x = rssringoccs_TwoVector_X(P);
    y = rssringoccs_TwoVector_Y(P);

    norm = rssringoccs_Sqrt_Double(x*x + y*y);
    return norm;
}

double rssringoccs_Euclidean_Norm_3D(rssringoccs_ThreeVector P)
{
    double x, y, z, norm;

    x = rssringoccs_ThreeVector_X(P);
    y = rssringoccs_ThreeVector_Y(P);
    z = rssringoccs_ThreeVector_Z(P);

    norm = rssringoccs_Sqrt_Double(x*x + y*y + z*z);
    return norm;
}

double
rssringoccs_Dot_Product_2D(rssringoccs_TwoVector P, rssringoccs_TwoVector Q)
{
    double Px, Py, Qx, Qy, dot_prod;

    Px = rssringoccs_TwoVector_X(P);
    Py = rssringoccs_TwoVector_Y(P);

    Qx = rssringoccs_TwoVector_X(Q);
    Qy = rssringoccs_TwoVector_Y(Q);

    dot_prod = Px*Qx + Py*Qy;
    return dot_prod;
}

double
rssringoccs_Dot_Product_3D(rssringoccs_ThreeVector P,
                           rssringoccs_ThreeVector Q)
{
    double Px, Py, Pz, Qx, Qy, Qz, dot_prod;

    Px = rssringoccs_ThreeVector_X(P);
    Py = rssringoccs_ThreeVector_Y(P);
    Pz = rssringoccs_ThreeVector_Z(P);

    Qx = rssringoccs_ThreeVector_X(Q);
    Qy = rssringoccs_ThreeVector_Y(Q);
    Qz = rssringoccs_ThreeVector_Z(Q);

    dot_prod = Px*Qx + Py*Qy + Pz*Qz;
    return dot_prod;
}

rssringoccs_TwoVector
rssringoccs_Normalize_TwoVector(rssringoccs_TwoVector P)
{
    double norm, rcpr_norm, x, y, x_hat, y_hat;
    rssringoccs_TwoVector P_normalized;

    x = rssringoccs_TwoVector_X(P);
    y = rssringoccs_TwoVector_Y(P);
    norm = rssringoccs_Euclidean_Norm_2D(P);
    rcpr_norm = 1.0/norm;

    x_hat = x*rcpr_norm;
    y_hat = y*rcpr_norm;

    P_normalized = rssringoccs_TwoVector_Rect(x_hat, y_hat);
    return P_normalized;
}

rssringoccs_ThreeVector
rssringoccs_Normalize_ThreeVector(rssringoccs_ThreeVector P)
{
    double norm, rcpr_norm, x, y, z, x_hat, y_hat, z_hat;
    rssringoccs_ThreeVector P_normalized;

    x = rssringoccs_ThreeVector_X(P);
    y = rssringoccs_ThreeVector_Y(P);
    z = rssringoccs_ThreeVector_Z(P);
    norm = rssringoccs_Euclidean_Norm_3D(P);
    rcpr_norm = 1.0/norm;

    x_hat = x*rcpr_norm;
    y_hat = y*rcpr_norm;
    z_hat = z*rcpr_norm;

    P_normalized = rssringoccs_ThreeVector_Rect(x_hat, y_hat, z_hat);
    return P_normalized;
}

rssringoccs_ThreeVector
rssringoccs_Cross_Product(rssringoccs_ThreeVector P, rssringoccs_ThreeVector Q)
{
    double Px, Py, Pz, Qx, Qy, Qz, x, y, z;
    rssringoccs_ThreeVector cross;

    Px = rssringoccs_ThreeVector_X(P);
    Py = rssringoccs_ThreeVector_Y(P);
    Pz = rssringoccs_ThreeVector_Z(P);

    Qx = rssringoccs_ThreeVector_X(Q);
    Qy = rssringoccs_ThreeVector_Y(Q);
    Qz = rssringoccs_ThreeVector_Z(Q);

    x = Py*Qz - Pz*Qy;
    y = Pz*Qx - Px*Qz;
    z = Px*Qy - Py*Qx;

    cross = rssringoccs_ThreeVector_Rect(x, y, z);
    return cross;
}



rssringoccs_TwoVector
rssringoccs_TwoVector_Add(rssringoccs_TwoVector P, rssringoccs_TwoVector Q)
{
    double Px, Py, Qx, Qy, x, y;
    rssringoccs_TwoVector sum;

    Px = rssringoccs_TwoVector_X(P);
    Py = rssringoccs_TwoVector_Y(P);

    Qx = rssringoccs_TwoVector_X(Q);
    Qy = rssringoccs_TwoVector_Y(Q);

    x = Px + Qx;
    y = Py + Qy;

    sum = rssringoccs_TwoVector_Rect(x, y);
    return sum;
}

rssringoccs_ThreeVector
rssringoccs_ThreeVector_Add(rssringoccs_ThreeVector P,
                            rssringoccs_ThreeVector Q)
{
    double Px, Py, Pz, Qx, Qy, Qz, x, y, z;
    rssringoccs_ThreeVector sum;

    Px = rssringoccs_ThreeVector_X(P);
    Py = rssringoccs_ThreeVector_Y(P);
    Pz = rssringoccs_ThreeVector_Z(P);

    Qx = rssringoccs_ThreeVector_X(Q);
    Qy = rssringoccs_ThreeVector_Y(Q);
    Qz = rssringoccs_ThreeVector_Z(Q);

    x = Px + Qx;
    y = Py + Qy;
    z = Pz + Qz;

    sum = rssringoccs_ThreeVector_Rect(x, y, z);
    return sum;
}

rssringoccs_TwoVector
rssringoccs_TwoVector_Scale(double a, rssringoccs_TwoVector P)
{
    double Px, Py, x, y;
    rssringoccs_TwoVector scale;

    Px = rssringoccs_TwoVector_X(P);
    Py = rssringoccs_TwoVector_Y(P);

    x = a*Px;
    y = a*Py;

    scale = rssringoccs_TwoVector_Rect(x, y);
    return scale;
}

rssringoccs_ThreeVector
rssringoccs_ThreeVector_Scale(double a, rssringoccs_ThreeVector P)
{
    double Px, Py, Pz, x, y, z;
    rssringoccs_ThreeVector scale;

    Px = rssringoccs_ThreeVector_X(P);
    Py = rssringoccs_ThreeVector_Y(P);
    Pz = rssringoccs_ThreeVector_Z(P);

    x = a*Px;
    y = a*Py;
    z = a*Pz;

    scale = rssringoccs_ThreeVector_Rect(x, y, z);
    return scale;
}

rssringoccs_ThreeVector
rssringoccs_Orthogonal_ThreeVector(rssringoccs_ThreeVector P)
{
    double Px, Py, x, y, z;
    rssringoccs_ThreeVector out;

    Px = rssringoccs_ThreeVector_X(P);
    Py = rssringoccs_ThreeVector_Y(P);

    if (Px == 0.0)
    {
        x = 1.0;
        y = 0.0;
        z = 0.0;
    }
    else
    {
        if (Py == 0.0)
        {
            x = 0.0;
            y = 1.0;
            z = 0.0;
        }
        else
        {
            x = 1.0;
            y = -Px/Py;
            z = 0.0;
        }
    }

    out = rssringoccs_ThreeVector_Rect(x, y, z);
    return out;
}
