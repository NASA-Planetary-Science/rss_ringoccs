"""

find_gaps.py

Purpose: Computes initial guess for radius of ring feature using the semimajor
         axis. Selects time and longitude closest to guess and computes true
         anomaly for a new radius guess. Continues this estimation method
         iteratively until difference between new and old radius guesses is
         less than some tolerance or maximum number of iterations is reached.

Notes:
    [1] Philip D. Nicholson, Richard G. French, Colleen A. McGhee-French, Matthew M. Hedman, Essam A. Marouf, Joshua E. Colwell, Katherine Lonergan, Talia Sepersky, Noncircular features in Saturn’s rings II: The C ring, Icarus, Volume 241, 2014, Pages 373-396, ISSN 0019-1035, https://doi.org/10.1016/j.icarus.2014.06.024. (http://www.sciencedirect.com/science/article/pii/S0019103514003443)
    [2] This must be run in a directory one level below the top-level rss_ringoccs directory, because of the csv file read in from tables/.
Revisions:
    2018 Oct 12 - sflury - original
    2017 Nov 19 - jfong - output both freespace_km and freespace_spm
"""

import numpy as np
import spiceypy as spice
import pdb

#def get_start_jd(year, doy):

def get_start_jd(year, doy):
    et = spice.utc2et(str(year)+'-'+(str(doy)).zfill(3)+'T00:00:00')
    jd = spice.et2utc(et, 'J', 10)
    start_jd = float(jd[3:])
    return start_jd


def find_gaps(t_ret_spm_vals, year, doy, rho_km_vals, phi_rl_deg_vals, niter=int(100),
        tolerance=0.001, t0=2454467.000000, kernels=None): #t0=252460865.18392503):

    if kernels:
        spice.kclear()
        spice.furnsh(kernels)

    # import tabulated info on gap geometry/orbits
    gaps = np.loadtxt('../tables/gap_orbital_elements.txt',
            skiprows=1, delimiter=',',usecols=(1,2,3,4))

    # places to store results
    gap_bounds = []

    t_start = get_start_jd(year, doy)
    
    time = t_start - t0 + (t_ret_spm_vals)/8.64e4
    

    # compute radii for all gaps
    for i in range(int(len(gaps)/2)):

        # storage for current gap's boundaries
        gap_limits = []

        # get inner and outer radii for current gap
        for j in range(2):

            # grab info for the gap's inner edge
            a,ae,w0,wdot = gaps[2*i+j]
            # compute excentricity
            e = ae/a
            # iteratively solve for radius
            radius = rad_converge( time, rho_km_vals, phi_rl_deg_vals,
                                a, e, w0, wdot, tolerance=tolerance,
                                niter=niter)
            # store in list
            gap_limits += [radius]

        # store limits for this gap in master list for return
        gap_bounds += [gap_limits]

    # return limits on all gaps
    return gap_bounds

# t0=2454467.000000 julian date time
# t0=252460865.18392503 ephemeris time
def rad_converge(t_oet_spm_vals,rho_km_vals,phi_rl_deg_vals,semimajor,
                eccentricity,curlypi_0,curlypi_dot, niter=int(100),
                tolerance=0.001):
    '''
    rad_converge

        Takes
            t_oet_spm_vals      -- ring intercept time (sec)
            rho_km_vals         -- ring intercept radius (km)
            phi_rl_deg_vals     -- interial longitude (deg/sec)
            semimajor           -- semimajor axis (km)
            eccentricity        -- eccentricity
            curlypi_0           -- longitude of periapse
            curlypi_dot         -- apsidal precession rate

        Optional
            t0        -- reference epoch (sec), default is UTC 12:00 Jan 1, 2008
            niter     -- maximum number of allowed iterations before terminating
            tolerance -- tolerance desired for solution convergence

        Computes initial guess for radius using the semimajor axis. Selects
            time and longitude closest to guess and computes true anomaly for
            a new radius guess. Continues this estimation method iteratively
            until difference between new and old radius guesses is less than
            some tolerance or the maximum number of iterations is reached.

        Returns
            radius_new      -- the approximate radius of the feature in question
                               computed by the iterative method
     '''
    # set initial guess at feature radius


    
    radius_old = semimajor*(1.+eccentricity)
    radius_guess = semimajor*(1.+eccentricity)

    # get index most closely associated with semi-major axis
    radius_diff = abs(rho_km_vals-radius_guess)
    ind = np.argwhere(radius_diff==min(radius_diff))
    # get associated time -- use median in case two values are equally good
    time = np.nanmedian(t_oet_spm_vals[ind])
    # get associated longitude -- used median in case two
    longitude = np.nanmedian(phi_rl_deg_vals[ind])
    # compute true anomaly from estimates of time and longitude using
    #    tabulated values of curly pi, curly pi dot, and reference epoch
    #    provided by the user
    true_anomaly = np.radians(longitude - curlypi_0 - curlypi_dot*(time))
    # use true anomaly to compute orbital radius from
    radius_new = (semimajor*(1-eccentricity**2))/(1+eccentricity*np.cos(true_anomaly))

    # iteration counter
    n = int(0)
    # continue computing new radius estimates until the difference between the
    #    previous and current (old and new) radius estimates are identical
    #    within some user-defined tolerance or the maximum number of iterations
    #    has been exceeded
    while abs(radius_new-radius_old) > tolerance and n < niter :

        # update radius estimate
        radius_guess = radius_new

        # get index most closely associated with semi-major axis
        radius_diff = abs(rho_km_vals-radius_guess)
        ind = np.argwhere(radius_diff==min(radius_diff))
        # get associated time -- use median in case two values are equally good
        time = np.nanmedian(t_oet_spm_vals[ind])
        # get associated longitude -- used median in case two
        longitude = np.nanmedian(phi_rl_deg_vals[ind])
        # compute true anomaly from estimates of time and longitude using
        #    tabulated values of curly pi, curly pi dot, and reference epoch
        #    provided by the user
        true_anomaly = np.radians(longitude - curlypi_0 - curlypi_dot*(time))
        # use true anomaly to compute orbital radius from
        radius_new = ((semimajor*(1-eccentricity**2))/
                     (1+eccentricity*np.cos(true_anomaly)))
        # update old radius estimate
        radius_old = radius_guess
        # update number of iterations
        n += 1

    return radius_new

