def diffraction_correction(data,res,out_file):
    import time as time
    t1 = time.time()
    import numpy as np
    from scipy.special import lambertw
    from scipy.special import iv
    import sys
    import spiceypy as spice
    import pandas as pd
    if type(data) != type('Hello'):
        sys.exit('ERROR: Input is not a string.')
    df = pd.read_csv(data,header=None,delim_whitespace=True)
    res_km           = np.double(res)
    rrange           = np.array([87410,87610])
    #Restore Variables
    rho_km_vals      = np.array(df[0])               #Ring intercept radius
    raw_power_vals   = np.array(df[2])               #Uncorrected power
    phase_rad_vals   = np.array(df[3])               #Uncorrected phase, rad
    b_rad_vals       = np.array(df[4])               #Ring-opening angle, rad
    d_km_vals        = np.array(df[5])               #Spacecraft-rip dist, km
    f_sky_hz_vals    = np.array(df[7])               #Sky frequency values, Hz
    phi_rad_vals     = np.array(df[8])               #Ring-azimuth angle, rad
    rho_dot_kms_vals = np.array(df[9])               #drho/dt, in km/s
    dx_km            = rho_km_vals[1]-rho_km_vals[0]
    drho_dt          = [np.min(rho_dot_kms_vals),np.max(rho_dot_kms_vals)]
    n                = len(rho_km_vals)               #Number of points
    n_iter           = 6                              #Number of iterations for psi
    dphi_rad         = 0.00001                        #Phase increment
     
    #If dx is negative, then occultation is ingress and variables must reverse.
    if dx_km < 0:
        rho_km_vals.reverse()
        phase_rad_vals.reverse()
        raw_power_vals.reverse()
        phi_rad_vals.reverse()
        b_rad_vals.reverse()
        f_sky_hz_vals.reverse()
        d_km_vals.reverse()
        rho_dot_kms_vals = -np.flipud(rho_dot_kms_vals) #dr/dt must be positive
        dx_km            = -dx_km
    
    #Compute other variables that are not contained in the original save file.
    lambda_sky  = spice.clight()/f_sky_hz_vals            #Wavelength of signal
    mu          = np.sin(np.abs(b_rad_vals))              #Sine of B
    f_km_vals   = np.sqrt(((lambda_sky*d_km_vals/2)*(1-(np.cos(b_rad_vals)**2)*\
    (np.sin(phi_rad_vals)**2)))/(np.sin(b_rad_vals)**2))  #The Fresnel Scale
    T_hat_vals  = np.sqrt(raw_power_vals)*np.exp(-1j*phase_rad_vals)
    sigma       = 2.e-13                                  #USO Allen Deviation
    omega       = 2.*np.pi*f_sky_hz_vals                  #Angular frequency (Hz)
    #Create variable to determine if an error has occurred.
    error_code = 0
    #Make sure the requested resolution has double valued precision.
    if (type(res_km) != type(0.)) and (type(res_km) != type(np.double(0))):
        print('Illegal res_km Type. Must be double valued.')
        error_code  +=  2
    
    #Check to make sure there is no ambiguity about the type of occultation.
    if (dx_km > 0) and (drho_dt[0] < 0):
        print("rho_km_vals are increasing but rho_dot_kms_vals is negative")
        print("I can't determine if this is an ingress or an egress occultation")
        print("Correct so rho_dot_kms_vals and rho_km_vals[1]-rho_km_vals[0]")
        print("Have the same sign.")
        error_code += 4
    
    #Check to make sure there is no ambiguity about the type of occultation.
    if (dx_km < 0) and (drho_dt[1] > 0):
        print("rho_km_vals are decreasing but rho_dot_kms_vals is positive")
        print("I can't determine if this is an ingress or an egress occultation")
        print("Correct so that rho_dot_kms_vals and rho_km_vals[1]-rho_km_vals[0]")
        print("Have the same sign.")
        error_code += 8
    
    #Check that the number of phase points equals the number of radius points.
    if len(phase_rad_vals) != n:
        print('Bad Input: len(phase_rad_vals) != len(rho_km_vals)')
        error_code +=  16
    
    #Check that the number of power points equals the number of radius points.
    if len(raw_power_vals) != n:
        print('Bad Input: len(raw_power_vals) != len(rho_km_vals)')
        error_code +=  32
    
    #Check the the ring azimuth and ring radius have the same number of points.
    if len(phi_rad_vals) != n:
        print('Bad Input: len(phi_rad_vals) != len(rho_km_vals)')
        error_code +=   64
    
    #Check that the ring-opening and ring radius have the same number of points.
    if len(b_rad_vals) != n:
        print('Bad Input: n_elements(B_rad_vals) ne n_elements(rho_km_vals)')
        error_code += 128
    
    #Check that the sky frequency and ring radius have the same number of points.
    if len(f_sky_hz_vals) != n:
        print('Bad Input: n_elements(f_sky_Hz_vals) ne n_elements(rho_km_vals)')
        error_code += 256
    
    #Check the the D-values and the ring radius have the same number of points.
    if len(d_km_vals) != n:
        print, 'Bad Input: n_elements(D_km_vals) ne n_elements(rho_km_vals)'
        error_code += 512
    
    #Check rho_dot_km_vals and the ring radius have the same number of points.
    if len(rho_dot_kms_vals) != n:
        print('Bad Input: n_elements(D_km_vals) ne n_elements(rho_km_vals)')
        error_code += 1024
    
    #Make sure that ring radius is a double precision array.
    if type(rho_km_vals) != type(np.zeros(10)):
        print('Bad Input: rho_km_vals is not double')
        error_code += 2048
    
    #Make sure that phase is a double precision array.
    if type(phase_rad_vals) != type(np.zeros(10)):
        print('Bad Input: phase_rad_vals is not double')
        error_code += 4096
    
    #Make sure that phase is a double precision array.
    if type(raw_power_vals) != type(np.zeros(10)):
        print('Bad Input: phase_rad_vals is not double')
        error_code += 8192
    
    #Make sure the ring-azimuth is a double precision array
    if type(phi_rad_vals) != type(np.zeros(10)):
        print('Bad Input: phi_rad_vals is not double')
        error_code += 16384
    
    #Make sure the B-values are a double precision array.
    if type(b_rad_vals) != type(np.zeros(10)):
        print('Bad Input: B_rad_vals is not double')
        error_code += 32678
    
    #Make sure the sky-frequency is a double precision array.
    if type(f_sky_hz_vals) != type(np.zeros(10)):
        print('Bad Input: f_sky_Hz_vals is not double')
        error_code += 65536
    
    #Make sure the D-values are a double precision array.
    if type(d_km_vals) != type(np.zeros(10)):
        print, 'Bad Input: D_km_vals is not double'
        error_code += 131072
    
    #Make sure the rho_dot_kms_vals are a double precision array.
    if type(rho_dot_kms_vals) != type(np.zeros(10)):
        print, 'Bad Input: rho_dot_kms_vals is not double'
        error_code += 262144
    
    #If a problem did occur, stop the code and find the error.
    if error_code != 0:
        sys.exit("Fatal Error Occured")
    
    #Parameters for inverse of Resolution = f(Window)
    alpha = (omega**2)*(sigma**2)/(2.*rho_dot_kms_vals)
    P     = res_km/(alpha*(f_km_vals**2))
    
    #The inverse exists only if P>1.
    if min(P) < 1.0001:
        print('WARNING: Bad Points!')
        print(np.where(P<1.0001))
        print('Either rho_dot_kms_vals, F_km_vals, or res_km is to small.')
        sys.exit("Fatal Error Occured")
    
    #Create window variable, window width (in km) for each point.
    Norm_eq   = 1.651920
    f         = np.abs(((P-1)*lambertw(np.exp(P/(1-P))*P/(1-P))+P)/(P-1))
    w_km_vals = Norm_eq*np.sqrt(2)*f/alpha
    w_max     = np.abs(np.max(w_km_vals))           #Largest window used, km
    
    #If the window width is large, inversion will take a while.
    if w_max > 1000:
        print('WARNING: Windows larger than 1000 km needed.')
        print('Max Window Width Used: ', w_max)
    
    #Smallest possible starting point.
    rho_min_lim   = np.min(rho_km_vals)+np.ceil(w_max/2.)
    #Largest possible ending point.
    rho_max_lim   = np.max(rho_km_vals)-np.ceil(w_max/2.)
    #Beggining of requested range
    rho_start     = rho_km_vals[np.min(np.where(rho_km_vals > rrange[0]))]
    #End of requested range
    rho_end       = rho_km_vals[np.max(np.where(rho_km_vals < rrange[1]))]
    #Starting radius, taking constraints into consideration.
    rho_min       = np.max([rho_min_lim,rho_start])
    #Ending radius, taking constraints into consideration.
    rho_max       = np.min([rho_max_lim,rho_end])
    #Starting point.
    start         = np.min(np.where(rho_km_vals > rho_min))
    #Ending point.
    finish        = np.max(np.where(rho_km_vals < rho_max))
    #Actual range used.
    range_km_vals = rho_km_vals[start:finish+1]
    #Number of points used.
    n_used        = len(range_km_vals)
    #Empty array for the corrected complex amplitude.
    T_vals        = np.zeros(n)*(1+1j)
    uni_m         = np.array([[1],[1],[1],[1],[1]])        #For matrix algebra
    #Create variable for Lagrange Interpolation
    dphi_val      = np.array([[-2.],[-1.],[0.],[1.],[2.]])*dphi_rad 
    #Calculate the corrected complex amplitude, point by point
    for i in range(0,n_used):
        center      = start+i
        n_pts_w     = int(2*np.ceil(w_km_vals[center]/(2*dx_km))-1)
        x           = np.array(range(0,n_pts_w))*dx_km - (n_pts_w-1)*dx_km/2.
        alph        = 2.5*np.pi
        w_func      = iv(0,alph*np.sqrt((1-(2*x/w_km_vals[center])**2)))/iv(0,alph)
        nw          = len(w_func)
        crange      = np.array(range(int(center-(nw-1)/2),int(1+center+(nw-1)/2)))
        r           = np.array([rho_km_vals[crange]])
        r0          = rho_km_vals[center]
        F_km_sq     = f_km_vals[center]**2
        D           = d_km_vals[center]
        B           = b_rad_vals[center]
        f_sky_hz    = f_sky_hz_vals[center]
        phi0        = phi_rad_vals[center]
        T_hat_Lu    = T_hat_vals[crange]
        lambda_km   = spice.clight()/f_sky_hz
        kD          = (spice.twopi()/lambda_km)*D
        cos_phi0    = np.cos(phi0)
        sin_phi0    = np.sin(phi0)
        sin_phi0_sq = sin_phi0**2
        r0cos_phi0  = r0*cos_phi0
        Dsq         = D**2
        r0_sq       = r0**2
        cosB        = np.cos(B)
        cosB_sq     = cosB**2
        factor      = (((cosB_sq)*cos_phi0*sin_phi0)/(1.-(cosB_sq)*sin_phi0_sq))/r0
        idd         = np.array([np.zeros(nw)+1])
        r_id        = np.dot(uni_m,r)
        r_sq        = np.dot(uni_m,r)**2
        phi_s_rad   = phi0-factor*(r-r0)
        dphi        = np.dot(dphi_val,idd)
        for kk in range(0,n_iter):
            phi_s_rad_vals = np.dot(uni_m,phi_s_rad) + dphi
            xi     = cosB*(r0cos_phi0-(r_id*np.cos(phi_s_rad_vals)))/D
            eta    = (r0_sq+r_sq-2.*r0*(r_id*np.cos(phi_s_rad_vals-phi0)))/Dsq
            FF     = kD*(np.sqrt(1.+2.*xi+eta)-(1.+xi))
            Psi_x  = (FF[0,::]-8.*FF[1,::]+8.*FF[3,::]-FF[4,::])
            Psi_xx = (-FF[0,::]+16.*FF[1,::]-30.*FF[2,::]+16.*FF[3,::]-FF[4,::])
            dphi_s_rad = -dphi_rad*Psi_x/Psi_xx
            phi_s_rad += dphi_s_rad
            
        xi              = (cosB/D)*(r0*cos_phi0-r*np.cos(phi_s_rad))
        eta             = (r0**2+r**2-2.*r*r0*np.cos(phi_s_rad-phi0))/Dsq
        psi_vals        = kD*(np.sqrt(1.+2.*xi+eta)-(1.+xi))
        psi_vals        = np.reshape(psi_vals,nw)
        arg             = np.exp(-1j*psi_vals)
        F_km            = f_km_vals[center]
        conv_kernel     = arg*w_func
        fft_t_hat       = np.fft.fft(T_hat_Lu)
        fft_conv        = np.fft.fft(conv_kernel)
        inv_t_hat       = np.fft.ifftshift(np.fft.ifft(fft_t_hat*fft_conv))
        inv_t_hat      *= dx_km*(np.complex(1.,1.))/(2.*F_km)
        T_vals[center]  = inv_t_hat[int((nw-1)/2)+1] 
        print(i, n_used - 1)
    
    r     = rho_km_vals[start:start+n_used]
    power = abs(T_vals[start:start+n_used])**2
    tau   = -mu[start:start+n_used]*np.log(power)
    t2    = time.time()
    np.savetxt(out_file,np.c_[r,power,tau],fmt='%32.16f '*3)
        
    print('Hours:', int(np.floor((t2-t1)/3600)))
    print('Minutes:',int(np.floor(t2-t1)/60))
    print('Seconds:', int(np.mod(t2-t1,60)))
