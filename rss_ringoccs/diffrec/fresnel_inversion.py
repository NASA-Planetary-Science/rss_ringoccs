def fresnel_transform_ellipse(T_in, rho_km_vals, F_km_vals, phi_rad_vals,
                              kD_vals, B_rad_vals, D_km_vals, w_km_vals, start,
                              n_used, peri, ecc, wtype, norm, fwd):

    # Compute the sample spacing.
    dx_km = rho_km_vals[1]-rho_km_vals[0]

    # Extract the window function from the window function dictionary.
    fw = window_functions.func_dict[wtype]["func"]

    # Create empty array for reconstruction / forward transform.
    T_out = T_in * 0.0

    # Compute first window width and window function.
    w_init = w_km_vals[start]

    # Number of points in the first window (Odd integer).
    nw = int(2 * np.floor(w_init / (2.0 * dx_km)) + 1)

    # Indices for the data corresponding to current window.
    crange = np.arange(int(start-(nw-1)/2), int(1+start+(nw-1)/2))

    # Various geometry variables.
    r0 = rho_km_vals[crange]
    r = rho_km_vals[start]
    x = r-r0

    # Compute the first window function.
    w_func = fw(x, w_init)

    # Perform the Fresnel Transform, point by point, via Riemann sums.
    for i in np.arange(n_used):

        # Current point being computed.
        center = start+i

        # Current window width, Fresnel scale, and ring radius.
        w = w_km_vals[center]
        F = F_km_vals[center]
        r = rho_km_vals[center]

        # If window widths changes too much, recompute the window function.
        if (np.abs(w_init - w) >= 2.0 * dx_km):

            # Compute first window width and window function.
            w_init = w_km_vals[center]

            # Number of points in window (Odd integer).
            nw = int(2 * np.floor(w_init / (2.0 * dx_km)) + 1)

            # Indices for the data corresponding to this window.
            crange = np.arange(int(center-(nw-1)/2), int(1+center+(nw-1)/2))

            # Ajdust ring radius by dx_km.
            r0 = rho_km_vals[crange]
            x = r-r0

            # Recompute the window function.
            w_func = fw(x, w_init)

        else:

            # If width hasn't changed much, increment index variable.
            crange += 1
            r0 = rho_km_vals[crange]

        # Various geometry variables for the current point.
        d = D_km_vals[crange]
        b = B_rad_vals[crange]
        kD = kD_vals[crange]
        phi = phi_rad_vals[crange]
        phi0 = phi_rad_vals[crange]

        # Compute Newton-Raphson perturbation
        psi_d1 = dpsi_ellipse(kD, r, r0, phi, phi0, b, d, ecc, peri)
        loop = 0
        while (np.max(np.abs(psi_d1)) > 1.0e-4):
            psi_d1 = dpsi_ellipse(kD, r, r0, phi, phi0, b, d, ecc, peri)
            psi_d2 = d2psi(kD, r, r0, phi, phi0, b, d)

            # Newton-Raphson
            phi += -(psi_d1 / psi_d2)

            # Add one to loop variable for each iteration
            loop += 1
            if (loop > 4):
                break

        # Compute Psi (Fresnel Kernel, MTR86 Equation 4).
        psi_vals = fresnel_psi(kD, r, r0, phi, phi0, b, d)

        # Compute kernel function for Fresnel inverse or forward model.
        if fwd:
            ker = w_func*np.exp(1j*psi_vals)
        else:
            ker = w_func*np.exp(-1j*psi_vals)

        # Compute approximate Fresnel transform for current point.
        T = T_in[crange]
        T_out[center] = np.sum(ker*T) * dx_km * (0.5+0.5j)/F

        # If normalization has been set, normalize the reconstruction
        if norm:
            T_out[center] *= window_functions.normalize(dx_km, ker, F)
    return T_out

def fresnel_transform_newton(T_in, rho_km_vals, F_km_vals, phi_rad_vals,
                             kD_vals, B_rad_vals, D_km_vals, w_km_vals, start,
                             n_used, wtype, norm, fwd):

    # Compute the sample spacing.
    dx_km = rho_km_vals[1]-rho_km_vals[0]

    # Extract the window function from the window function dictionary.
    fw = window_functions.func_dict[wtype]["func"]

    # Create empty array for reconstruction / forward transform.
    T_out = T_in * 0.0

    # Compute first window width and window function.
    w_init = w_km_vals[start]

    # Number of points in the first window (Odd integer).
    nw = int(2 * np.floor(w_init / (2.0 * dx_km)) + 1)

    # Indices for the data corresponding to current window.
    crange = np.arange(int(start-(nw-1)/2), int(1+start+(nw-1)/2))

    # Various geometry variables.
    r0 = rho_km_vals[crange]
    r = rho_km_vals[start]
    x = r-r0

    # Compute the first window function.
    w_func = fw(x, w_init)

    # Perform the Fresnel Transform, point by point, via Riemann sums.
    for i in np.arange(n_used):

        # Current point being computed.
        center = start+i

        # Current window width, Fresnel scale, and ring radius.
        w = w_km_vals[center]
        F = F_km_vals[center]
        r = rho_km_vals[center]

        # If window widths changes too much, recompute the window function.
        if (np.abs(w_init - w) >= 2.0 * dx_km):

            # Compute first window width and window function.
            w_init = w_km_vals[center]

            # Number of points in window (Odd integer).
            nw = int(2 * np.floor(w_init / (2.0 * dx_km)) + 1)

            # Indices for the data corresponding to this window.
            crange = np.arange(int(center-(nw-1)/2), int(1+center+(nw-1)/2))

            # Ajdust ring radius by dx_km.
            r0 = rho_km_vals[crange]
            x = r-r0

            # Recompute the window function.
            w_func = fw(x, w_init)
        else:

            # If width hasn't changed much, increment index variable.
            crange += 1
            r0 = rho_km_vals[crange]

        # Various geometry variables for the current point.
        d = D_km_vals[crange]
        b = B_rad_vals[crange]
        kD = kD_vals[crange]
        phi = phi_rad_vals[crange]
        phi0 = phi_rad_vals[crange]

        # Compute dpsi/dphi.
        psi_d1 = fresnel_dpsi_dphi(kD, r, r0, phi, phi0, b, d)

        # Variable for breaking out of the inner for loop.
        loop = 0

        # Perform Newton-Raphson to find the roots of dpsi/dphi.
        while (np.max(np.abs(psi_d1)) > 1.0e-4):

            # Compute the first and second derivatives of psi.
            psi_d1 = fresnel_dpsi_dphi(kD, r, r0, phi, phi0, b, d)
            psi_d2 = d2psi(kD, r, r0, phi, phi0, b, d)

            # Newton-Raphson Perturbation.
            phi += -(psi_d1 / psi_d2)

            # Add one to loop variable for each iteration
            loop += 1
            if (loop > 4):
                break

        # Compute Psi (Fresnel Kernel, MTR86 Equation 4).
        psi_vals = fresnel_psi(kD, r, r0, phi, phi0, b, d)

        # Compute kernel function for Fresnel inverse or forward model.
        if fwd:
            ker = w_func*np.exp(1j*psi_vals)
        else:
            ker = w_func*np.exp(-1j*psi_vals)

        # Compute approximate Fresnel transform for current point.
        T = T_in[crange]
        T_out[center] = np.sum(ker*T) * dx_km * (0.5+0.5j)/F

        # If normalization has been set, normalize the reconstruction
        if norm:
            T_out[center] *= window_functions.normalize(dx_km, ker, F)
    return T_out

def fresnel_transform_quadratic(T_in, rho_km_vals, F_km_vals, w_km_vals, start,
                                n_used, wtype, norm, fwd):

    # Compute the sample spacing.
    dx_km = rho_km_vals[1]-rho_km_vals[0]

    # Extract the window function from the window dictionary.
    fw = window_functions.func_dict[wtype]["func"]

    # Create empty array for reconstruction / forward transform.
    T_out = T_in * 0.0

    # Extract the first window width.
    w_init = w_km_vals[start]

    # Compute the number of points needed for window (Must be odd integer).
    nw = int(2 * np.floor(w_init / (2.0 * dx_km)) + 1)

    # Range of indices of rho_km_vals that are processed for first window.
    crange = np.arange(int(start-(nw-1)/2), int(1+start+(nw-1)/2))

    # Geometry and data to reconstruct the first point.
    r0 = rho_km_vals[crange]
    r = rho_km_vals[start]
    x = r-r0
    x2 = HALF_PI * np.square(x)

    # Compute the window function about the current point.
    w_func = fw(x, w_init)

    # Decrement crange by 1 since a +1 occurs in the next for loop.
    crange -= 1

    # Precompute the square of F_km_vals for speed.
    F2 = np.square(F_km_vals)

    # Perform Fresnel Transform, point by point, using Riemann Sums.
    for i in np.arange(n_used):

        # Current point being computed.
        center = start+i

        # Window width and Frensel scale for current point.
        w = w_km_vals[center]
        F = F_km_vals[center]

        # If the width has changed too much, recompute the window function.
        if (np.abs(w_init - w) >= 2.0 * dx_km):

            # Compute first window width and window function.
            w_init = w_km_vals[center]

            # Number of points in window (Odd integer).
            nw = int(2 * np.floor(w_init / (2.0 * dx_km))+1)

            # Range of indices for rho_km_vals for this window.
            crange = np.arange(int(center-(nw-1)/2), int(1+center+(nw-1)/2))

            # Ajdust ring radius by dx_km.
            r0 = rho_km_vals[crange]
            r = rho_km_vals[center]

            # Center rho_km_vals about r0.
            x = r-r0
            x2 = HALF_PI * np.square(x)

            # Recompute the window function.
            w_func = fw(x, w_init)
        else:

            # If width hasn't changed much, increment indexing variable.
            crange += 1

        # Fresnel approximation gives psi = pi/2 ((r-r0) / F)^2.
        psi_vals = x2 / F2[center]

        # Compute kernel function for Fresnel inverse or forward model.
        if fwd:
            ker = w_func*np.exp(1j*psi_vals)
        else:
            ker = w_func*np.exp(-1j*psi_vals)

        # Compute approximate Fresnel transform for current point
        T_out[center] = np.sum(ker * T_in[crange]) * dx_km * (0.5+0.5j) / F

        # If normalization has been set, normalize the reconstruction.
        if norm:
            T_out[center] *= window_functions.normalize(dx_km, ker, F)
    return T_out

def fresnel_legendre_transform(T_in, rho_km_vals, F_km_vals, phi_rad_vals,
                               kD_vals, B_rad_vals, D_km_vals, w_km_vals, start,
                               n_used, wtype, norm, fwd, psitype):

    if (psitype == "fresnel4"):
        ord = 3
    elif (psitype == "fresnel6"):
        ord = 5
    elif (psitype == "fresnel8"):
        ord = 7
    else:
        raise ValueError(
            "\r\tError Encountered: rss_ringoccs"
            "\r\t\tdiffrec.special_functions.fresnel_legendre_transform\n"
            "\r\tpsitype must be set to 'fresnel4', 'fresnel6', or 'fresnel8'"
        )

    dx_km = rho_km_vals[1]-rho_km_vals[0]

    # Define functions.
    fw = window_functions.func_dict[wtype]["func"]

    # Create empty array for reconstruction / forward transform.
    T_out = T_in * 0.0

    # Compute first window width.
    w_init = w_km_vals[start]

    # Compute the number of points needed for window (Must be odd integer).
    nw = int(2 * np.floor(w_init / (2.0 * dx_km)) + 1)

    # Range of indices of rho_km_vals that are processed for first window.
    crange = np.arange(int(start-(nw-1)/2), int(1+start+(nw-1)/2))

    # Geometry and data to reconstruct the first point.
    r0 = rho_km_vals[crange]
    r = rho_km_vals[start]
    x = r-r0
    x2 = np.square(x)

    # Window function about first point.
    w_func = fw(x, w_init)

    # Decrement crange by 1 since a +1 occurs in the next for loop.
    crange -= 1

    # Precompute sine and cosine of variables for speed.
    cosb = np.cos(B_rad_vals)
    cosp = np.cos(phi_rad_vals)
    sinp = np.sin(phi_rad_vals)
    Legendre_Coeff  = cosb*sinp
    Legendre_Coeff *= Legendre_Coeff
    Legendre_Coeff  = 0.5*Legendre_Coeff/(1.0-Legendre_Coeff)

    # Compute coefficient that occurs in Legendre expansion of psi.
    A_2 = 0.5*np.square(cosb*sinp)/(1.0-np.square(cosb*sinp))

    # Legendre Polynomials.
    legendre_p = []
    legendre_p.append(1.0)
    legendre_p.append(cosb*cosp)
    fresnel_p = []
    fresnel_p.append(0.5-0.5*np.square(legendre_p[1]))
    coeffs_p = []
    l_coeffs = []
    for i in range(ord):
        l_coeffs.append(1.0/(i+2.0))

    for i in range(1, ord):
        legendre_p.append(((2.0*i+1.0)*legendre_p[1]*legendre_p[i] -
                          i*legendre_p[i-1])*l_coeffs[i-1])
        fresnel_p.append((legendre_p[i]-legendre_p[1]*legendre_p[i+1])*l_coeffs[i])

    for i in range(1, int((ord+1)/2)+1):
        coeffs_p.append(0.0)
        for k in range(i):
            coeffs_p[i-1] += legendre_p[k+1]*legendre_p[i-k]
        coeffs_p[i-1] = fresnel_p[i-1] - Legendre_Coeff*coeffs_p[i-1]

    for i in range(int((ord+1)/2)+1, ord):
        coeffs_p.append(0.0)
        for k in range(i - int((ord+1)/2), int((ord+1)/2)):
            coeffs_p[i-1] += legendre_p[k+1]*legendre_p[i-k]
        coeffs_p[i-1] = fresnel_p[i-1] - Legendre_Coeff*coeffs_p[i-1]

    i = int((ord+1)/2)
    coeffs_p.append((fresnel_p[ord-1]-Legendre_Coeff*np.square(legendre_p[i])))
    rcpr_d = 1.0/D_km_vals
    rcpr_d2 = np.square(rcpr_d)

    # Perform Fresnel transform, point by point, using Riemann sums.
    for i in np.arange(n_used):

        # Current point being computed.
        center = start+i

        # Window width and Frensel scale for current point.
        w = w_km_vals[center]
        F = F_km_vals[center]

        # If the width has changed too much, recompute the window function.
        if (np.abs(w_init - w) >= 2.0 * dx_km):

            # Compute first window width and window function.
            w_init = w_km_vals[center]

            # Number of points in window (Odd integer).
            nw = int(2 * np.floor(w_init / (2.0 * dx_km))+1)

            # Range of indices for rho_km_vals for this window.
            crange = np.arange(int(center-(nw-1)/2), int(1+center+(nw-1)/2))

            # Ajdust ring radius by dx_km.
            r0 = rho_km_vals[crange]
            r = rho_km_vals[center]

            # Center rho_km_vals about r0.
            x = r-r0
            x2 = np.square(x)

            # Recompute the window function.
            w_func = fw(x, w_init)
        else:

            # If width hasn't changed much, increment indexing variable.
            crange += 1

        # Independent variable used in the Legendre expansion.
        z = x*rcpr_d[center]
        psi_vals = coeffs_p[ord-1][center]
        for k in range(2, ord):
            psi_vals = psi_vals*z + coeffs_p[ord-k][center]
        
        psi_vals = psi_vals*z + coeffs_p[0][center]
        psi_vals *= kD_vals[center] * np.square(z)

        # Compute kernel function for Fresnel inverse or forward model.
        if fwd:
            ker = w_func*np.exp(1j*psi_vals)
        else:
            ker = w_func*np.exp(-1j*psi_vals)

        # Compute approximate Fresnel transform for current point
        T_out[center] = np.sum(ker * T_in[crange]) * dx_km * (0.5+0.5j) / F

        # If normalization has been set, normalize the reconstruction.
        if norm:
            T_out[center] *= window_functions.normalize(dx_km, ker, F)
    return T_out
