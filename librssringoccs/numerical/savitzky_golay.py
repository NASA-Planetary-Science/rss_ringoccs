/*  You'll need these codes later. */
static void _Moore_Penrose_Pseudo_Inverse;
static void _Singular_Value_Decomposition;

_Bool Savitzky_Golay_Double(double *x_in, double *x_smooth, long arr_size,
                            long window_size, int deriv_order, int poly_order)
{
    if (((window_size % 2) != 1) || (window_size < 1))
        return 0;
    else if (window_size < poly_order + 2)
        return 0;

    long half_window = (window_size - 1)/2;

    int n, m, pos_K, minus_K, K;
    long factor;

    /*  Precompute the Savitzky-Golay coefficients.                           */
    double **mat;

    mat = malloc(sizeof(*b) * window_size);

    for (n=0; n<window_size; ++n)
        mat[n] = calloc(sizeof(*mat[n]), poly_order);

    for (n=0; n<window_size; ++n)
        mat[n][0] = 1.0;

    for (n=0; n<half_window; ++n)
    {
        K       = half_window - n;
        pos_K   = K;
        minus_K = -K;

        for (m=1; m<=poly_order; ++m)
        {
            mat[n][m] = minus_K;
            mat[window_size-1-n][m] = pos_K;
            minus_K *= -K;
            pos_K   *= K;
        }
    }

    for (m=0; m<deriv_order; ++m)
        factor *= (deriv_order-m);

    for k in range(half_window):
        n0 = half_window - k
        m = n0
        n = -n0
        for j in range(1, order+1):
            b[k, j] = n
            b[window_size-1-k, j] = m
            n *= -n0
            m *= n0

    b = numpy.mat(b)

    m = 1
    for i in range(deriv):
        m *= rate*(deriv-i)

    m *= numpy.linalg.pinv(b).A[deriv]

    # Pad the endpoints with values from the signal.
    firstvals = y[0] - numpy.abs(y[1:half_window+1][::-1] - y[0])
    lastvals = y[-1] + numpy.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = numpy.concatenate((firstvals, y, lastvals))

    return numpy.convolve(m[::-1], y, mode='valid')
}