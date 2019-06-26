"""

Purpose: Create the plots in the second version of Essam's EASYDATA summary PDF.

"""

import warnings
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib
matplotlib.rcParams['agg.path.chunksize'] = 10000
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.signal import savgol_filter
from .sys_tools import latex_summary_doc
from .write_output_files import construct_filepath
import time
import spiceypy as spice

sat_radius =  60268.
rings_km = [74490., 91983., 117516., 122052., 136774., 139826.]
radii = [74.490, 91.983, 117.516, 122.052, 133.424, 136.774, 140.461]
Saturn = 699
Earth   = 399
Cassini = -82
ref = 'J2000'
planet = Saturn


Rvals=[74490., 91983.,117516., 122052.,136774.]
lw1 = 1.0

def plot_bullseye(pdf, dlp_inst):
    # Grab geometry information
    plt.close()
    rho_km = dlp_inst.rho_km_vals
    phi_rl_rad = dlp_inst.phi_rl_rad_vals
    phi_ora_rad = dlp_inst.phi_rad_vals
    B_mean = round(np.mean(np.degrees(dlp_inst.B_rad_vals)), 4)
    band = dlp_inst.rev_info['band']
    t_oet_spm = dlp_inst.t_oet_spm_vals


    # Radio track relative to Earth direction
    occ_ed_x = rho_km * np.cos(phi_ora_rad)
    occ_ed_y = rho_km * np.sin(phi_ora_rad)
  
    # Radio track relative to inertial referenc
    occ_an_x = rho_km * np.cos(phi_rl_rad)
    occ_an_y = rho_km * np.sin(phi_rl_rad)
    
    theta = np.linspace(0., 2.*np.pi, 361)
    plt.figure()

    xx = np.cos(theta)
    yy = np.sin(theta)

    # Plot Saturn
    sat_xx = sat_radius/1000. * xx
    sat_yy = sat_radius/1000. * yy
    plt.plot(sat_xx, sat_yy, 'k-', linewidth=1.5)
    plt.plot(0,0,'+')
    plt.text(5,-3, 'Saturn')

    # Plot occultation tracks
    plt.plot(occ_an_x/1000., occ_an_y/1000., 'b-', linewidth=1.5)
    plt.plot(occ_ed_x/1000., occ_ed_y/1000., 'r-', linewidth=1.5)
    plt.axis('equal')
    plt.xlabel('$x$ (1000 km)')
    plt.ylabel('$y$ (1000 km)')
    plt.title('Mean Observed Ring Elevation $B$ = '+str(B_mean)+'$^\circ$')
    plt.tick_params(axis='both', direction='in', top=True, bottom=True, right=True, left=True)

    # Plot rings
    for ring in rings_km:
        plt.plot(ring/1000.*xx, ring/1000.*yy, 'm-', linewidth=1.)
    plt.text(80, -3, 'C')
    plt.text(101, -3, 'B')
    plt.text(125, -3, 'A')
    plt.xlim([-150., 150.])
    plt.ylim([-150., 150.])

    pdf.savefig()
    plt.close()
    return pdf
def calc_bp(rap,decp,ra,dec):
    # Compute B and P from Eqs A49-51 of 28 Sgr geometry paper
    dtor=np.pi/180.
    ap=rap*dtor
    dp=decp*dtor
    a=ra*dtor
    d=dec*dtor
    Brad=np.arcsin(-np.cos(dp)*np.cos(d)*np.cos(a-ap) - np.sin(dp)*np.sin(d))
    B=Brad/dtor
    Prad=np.arctan([-np.cos(dp)*np.sin(a-ap)/(-np.cos(dp)*np.sin(d)*np.cos(a-ap)+np.sin(dp)*np.cos(d))])
    P=Prad/dtor
    return B, P

# note: p2 is a legend page with no plots
def plot_occ_earth_view(pdf, geo_inst):

    # initialize plot
    plt.close()
    fig = plt.figure(1)
    plt.axis('off')
    plt.axes().set_aspect('equal')

    kernels = geo_inst.history['Positional Args']['kernels']
    spice.kclear()
    spice.furnsh(kernels)

    t_oet_et_vals = geo_inst.t_oet_et_vals
    oet_et_start = t_oet_et_vals[0]
    t_set_et_vals = geo_inst.t_set_et_vals
    rho_km = geo_inst.rho_km_vals
    phi_rl_deg = geo_inst.phi_ora_deg_vals
    set_et_start = t_set_et_vals[0]
    # Compute transformation between J2000 and planet coordinates
    #   EEQ = Earth Equatorial System
    #   PEQ = Planetary ephemeris system
    mm_EEQ2PEQ = spice.tipbod(ref, planet, set_et_start)
    mm_PEQ2EEQ = np.transpose(mm_EEQ2PEQ)

    # planet pole coordinates in J2000
    nhat_pole_PEQ = [0., 0., 1.]


    nhat_pole_planet = np.dot(mm_PEQ2EEQ, nhat_pole_PEQ)
    radius, rapole_planet, decpole_planet = spice.recrad(nhat_pole_planet)


    # compute earth to planet vector at start of occ
    ptarg, ltime = spice.spkezp(planet, oet_et_start, ref, 'LT', Earth)
    radius_planet, ra_planet, dec_planet = spice.recrad(ptarg)

    rap = rapole_planet*spice.dpr() # in degrees
    decp= decpole_planet*spice.dpr() # in degrees
    ra  = ra_planet*spice.dpr() # in degrees
    dec = dec_planet*spice.dpr() # in degrees

    Bdeg, Pdeg = calc_bp(rap, decp, ra, dec)

    # transformations between EEQ, PEQ and sky plane
    Prad = Pdeg * spice.rpd()
    mm_SKY2EEQ = spice.eul2m(-ra_planet, dec_planet, Prad, 3, 2, 1)
    mm_EEQ2SKY = np.transpose(mm_SKY2EEQ)

    mm_PEQ2SKY = np.matmul(mm_EEQ2SKY, mm_PEQ2EEQ)
    mm_SKY2PEQ = np.transpose(mm_PEQ2SKY)

    # get equatorial and polar radii of planet
    abc = spice.bodvar(planet, 'RADII', 3)
    Re = abc[0]
    Rp = abc[2]

    # compute projected polar radius in sky plane
    Brad = Bdeg * spice.rpd()
    Rpsky = np.sqrt((Re*np.sin(Brad))**2 + (Rp*np.cos(Brad))**2)

    # draw outline of planet on sky, rotated by position angle of pole
    ntheta = 361
    theta = np.linspace(0., 360., ntheta)
    cosTheta = np.cos(theta*spice.rpd())
    sinTheta = np.sin(theta*spice.rpd())

    u = Re * cosTheta
    v = Rpsky * sinTheta
    cosP = np.cos(Prad) * 0+1.
    sinP = np.sin(Prad) *0
    x = cosP * u - sinP * v
    y = sinP * u + cosP * v

    plt.plot(x, y, color='black', linewidth=lw1)

    # Add rings
    xsky = np.zeros(ntheta)
    ysky = np.zeros(ntheta)
    zsky = np.zeros(ntheta)
    for nr in range(len(Rvals)):
        R = Rvals[nr]
        Rx = R * cosTheta
        Ry = R * sinTheta
        Rz = 0.
        for n in range(ntheta):
            vec1 = [Rx[n], Ry[n], Rz]
            vecsky = np.dot(mm_PEQ2SKY, vec1)
            xsky[n] = vecsky[1]
            ysky[n] = vecsky[2]
            zsky[n] = vecsky[0]

        Rsky = np.sqrt(xsky**2 + (Re*ysky/Rpsky)**2)
        req1 = zsky < 0
        req2 = Rsky > Re
        l = np.logical_or(req1, req2)
        nl = np.sum(l)
        xsky[~l] = np.nan
        ysky[~l] = np.nan
        plt.plot(-xsky, ysky, color='red', linewidth=lw1)

    # Add earth direction
    plt.plot([0.,0.], [-6.1e4*np.sin(Brad), -1.8e5*np.sin(Brad)],
            color='black', linewidth=2.0)
    plt.text(500, -1.6e5*np.sin(Brad), '$\oplus$', fontsize=15)


    # add latitude circles every 10 degrees
    nlats = 16.
    latitudes = np.radians((-80.0 + np.arange(nlats) * 10.))  # in radians
    for nlat in range(int(nlats)):
        lat = latitudes[nlat]
        R1 = Re * np.cos(lat)
        Rx = R1 * cosTheta
        Ry = R1 * sinTheta
        Rz = Rpsky * np.sin(lat)
        for n in range(ntheta):
            vec1 = [Rx[n], Ry[n], Rz]
            vecsky = np.dot(mm_PEQ2SKY, vec1)
            xsky[n] = vecsky[1]
            ysky[n] = vecsky[2]
            zsky[n] = vecsky[0]

        Rsky = np.sqrt(xsky**2 + (Re*ysky/Rpsky)**2)
        l = zsky < 0
        xsky[~l] = np.nan
        ysky[~l] = np.nan
        plt.plot(-xsky, ysky, color='black', linewidth=lw1)



    # add meridians every 20 degrees
    nphi = 181.
    phi = np.radians(np.arange(nphi))       # in radians
    nlons = 18.
    longitudes = np.radians((-180. + np.arange(nlons) * 20.))      # in radians
    zsat = Rpsky * np.cos(phi)
    nphi_plot = []

    # loop over meridians
    for nlon in range(int(nlons)):
        longitude = longitudes[nlon]
        xsat = Re * np.cos(longitude) * np.sin(phi)
        ysat = Re * np.sin(longitude) * np.sin(phi)
        nplot = []
        for n in range(1,int(nphi)):
            vecsat = [xsat[n], ysat[n], zsat[n]]
            vecsky = np.dot(mm_PEQ2SKY, vecsat)
            xskyval = vecsky[1]
            yskyval = vecsky[2]
            zskyval = vecsky[0]
            if zskyval < 0.:
                xvalp = -xskyval
                yvalp = yskyval
                nplot.append([xvalp, yvalp])

        nphi_plot.append(nplot)

    for plist in nphi_plot:
        xarr = []
        yarr = []
        for arr in plist:
            xarr.append(arr[0])
            yarr.append(arr[1])
        plt.plot(xarr, yarr, color='black', linewidth=lw1)


    # add occultation track
    rho_mask = np.array([True for x in range(len(rho_km))])

    for ind in range(len(rho_km)):
        if rho_km[ind] > 150000.:
            rho_mask[ind] = False
    occx1 = rho_km[rho_mask] * np.cos(np.radians(phi_rl_deg[rho_mask]-90.)) * np.cos(Brad)
    occy1 = rho_km[rho_mask]* np.sin(np.radians(phi_rl_deg[rho_mask]-90.)) * np.sin(Brad)
    blocked_mask = np.array([False for x in range(len(occx1))])
    for i in range(len(occx1)):
        if geo_inst.t_oet_et_vals[i] in geo_inst.ionos_occ_et_vals:
            blocked_mask[i] = True
    plt.plot(occx1[blocked_mask], occy1[blocked_mask], linestyle='--',
            color='blue')
    occx1[blocked_mask] = np.nan
    occy1[blocked_mask] = np.nan
    plt.plot(occx1, occy1, color='blue')
    dpts = 60 * 30
    plt.scatter(occx1[::dpts], occy1[::dpts], s=10, color='blue')
        
    pdf.savefig()
    plt.close()

    return pdf

def plot_occ_pole_view(pdf, geo_inst):
    plt.close()
    plt.axis('off')
    plt.axes().set_aspect('equal')

    phi_ora_deg = geo_inst.phi_ora_deg_vals
    rho_km = geo_inst.rho_km_vals

    # draw outline of planet 
    ntheta = 361
    theta = np.linspace(0., 360., ntheta)
    cosTheta = np.cos(np.radians(theta))
    sinTheta = np.sin(np.radians(theta))
    abc = spice.bodvar(planet, 'RADII', 3)
    Re = abc[0]
    Rp = abc[1]

    u = Re * cosTheta
    v = Rp * sinTheta
    cosP = np.cos(np.pi/2) * 0+1.
    sinP = np.sin(np.pi/2) *0
    x = cosP * u - sinP * v
    y = sinP * u + cosP * v

    plt.plot(x, y, color='black', linewidth=lw1)

    # add rings
    for nr in range(len(Rvals)):
        R = Rvals[nr]
        Rx = R * cosTheta
        Ry = R * sinTheta
        plt.plot(Rx, Ry, color='red', linewidth=0.5)

    # add latitude markers
    nlats = 9
    latitudes = np.radians(np.arange(nlats)*10.)
    #latitudes = np.radians((-80.0 + np.arange(nlats) * 10.))  # in radians
    for nlat in range(int(nlats)):
        lat = latitudes[nlat]
        R1 = Re * np.cos(lat)
        Rx = R1 * cosTheta
        Ry = R1 * sinTheta
        Rz = Rp * np.sin(lat)
        plt.plot(Rx, Ry, color='black', linewidth=lw1)

    # add longitude markers
    nlongs = 18
    longitudes = np.radians(np.arange(nlongs)*20.)
    for nlong in longitudes:
        rx = Re * np.cos(nlong)
        ry = Re * np.sin(nlong)
        plt.plot([0.,rx], [0.,ry], color='black', linewidth=lw1)

    rho_mask = np.array([True for x in range(len(rho_km))])

    for ind in range(len(rho_km)):
        if rho_km[ind] > 150000.:
            rho_mask[ind] = False

    # add occultation track
    occx = rho_km[rho_mask] * np.cos(np.radians(phi_ora_deg[rho_mask]-90.))
    occy = rho_km[rho_mask]* np.sin(np.radians(phi_ora_deg[rho_mask]-90.))

    occr = np.sqrt(occx**2 + occy**2)
    plt.plot(occx[occr < Re], occy[occr<Re], linestyle='--', color='blue',
             linewidth=lw1)

    occx = occx[occr > Re]
    occy = occy[occr > Re]


    plt.plot(occx, occy, color='blue')

    # add 30-min markers, assuming each point is spaced 1sec apart
    dpts = 60 * 30
    plt.scatter(occx[::dpts], occy[::dpts], color='blue', s=10.)

    # add direction to earth
    plt.plot([0.,0.], [-68000., -180000.], 'k-', linewidth=2.3)
    plt.text(2000., -180000., '$\oplus$', fontsize=15)
    pdf.savefig()
    plt.close()
    return pdf

def plot_geo_overview(pdf, geo_inst, tau_inst):
    t_oet_1kspm = geo_inst.t_oet_spm_vals/1000.
    t_set_1kspm = geo_inst.t_set_spm_vals/1000.
    rho_km = geo_inst.rho_km_vals
    phi_rl_deg = geo_inst.phi_rl_deg_vals
    phi_ora_deg = geo_inst.phi_ora_deg_vals
    rho_dot_kms = geo_inst.rho_dot_kms_vals
    D_km = geo_inst.D_km_vals
    F_km = geo_inst.F_km_vals
    band = geo_inst.rev_info['band']
    phi_rl_dot_kms = geo_inst.phi_rl_dot_kms_vals
    mask_rad = np.array([True for x in range(len(rho_km))])
    mask_az= np.array([True for x in range(len(rho_km))])
    for val in range(len(rho_km)):
        if rho_dot_kms[val] < -150 or rho_dot_kms[val] > 150:
            mask_rad[val] = False
        if phi_rl_dot_kms[val] < -150 or phi_rl_dot_kms[val] > 150:
            mask_az[val] = False



    tau_threshold_vals = tau_inst.tau_threshold_vals
    rho_km_tau = tau_inst.rho_km_vals

    nrow = 3
    ncol = 2
    fig, axes = plt.subplots(nrow, ncol, figsize=(8.5, 11))
    ax1 = axes[0,0]
    ax2 = axes[0,1]
    ax3 = axes[1,0]
    ax4 = axes[1,1]
    ax5 = axes[2,0]
    ax6 = axes[2,1]

    # gridspec_kw = {'width_ratios':[3, 1]}
    
    
    ax1.set_title('Time (1000 SPM ' + geo_inst.rev_info['year'] + '-'
                + geo_inst.rev_info['doy'] +')', fontweight='bold', fontsize=10)
    ytpos = 0.08
    ax1.plot(rho_km/1000., t_oet_1kspm, 'b', label = 'ERT')
    ax1.plot(rho_km/1000., t_set_1kspm, 'r', label='SCET')
    ax1.text(0.17, ytpos, 'C', horizontalalignment='center',
            verticalalignment='center', transform=ax1.transAxes)

    ax1.text(0.43, ytpos, 'B', horizontalalignment='center',
            verticalalignment='center', transform=ax1.transAxes)
    ax1.text(0.75, ytpos, 'A', horizontalalignment='center',
            verticalalignment='center', transform=ax1.transAxes)
    ax1.text(0.94, ytpos, 'F', horizontalalignment='center',
            verticalalignment='center', transform=ax1.transAxes)
    ax1.legend(loc='upper right')

    
    ax2.set_title('S/C Distance D/RS', fontweight='bold', fontsize=10)
    ax2.plot(rho_km/1000., D_km/sat_radius, 'k')
    
    ax3.set_title('Longitude, Azimuth (deg)',fontweight='bold', fontsize=10)
    ax3.plot(rho_km/1000., phi_rl_deg, 'b', label = '$\phi_{J2K}$')
    ax3.plot(rho_km/1000., phi_ora_deg, 'r', label = '$\phi_{E}$')
    ax3.legend(loc='upper right')
    
    
    ax4.set_title('Occ Track Velocity (km/s)', fontweight='bold', 
             fontsize=10)
    ax4.plot(rho_km[mask_rad]/1000., rho_dot_kms[mask_rad], 'r', label = '$V_{rad}$')
    ax4.plot(rho_km[mask_az]/1000., phi_rl_dot_kms[mask_az], 'b', label = '$V_{az}$')
    ax4.legend(loc='upper right')
    
    ax5.set_title('Fresnel Scale (km)', fontweight='bold', fontsize=10)
    ax5.plot(rho_km/1000., F_km, 'r', label = band.split('"')[1])
    ax5.legend(loc='upper right')
    
    ax6.set_title('Threshold Optical Depth', fontweight='bold',
            fontsize=10)
    ax6.plot(rho_km_tau/1000., tau_threshold_vals, label=band.split('"')[1])
    ax6.legend(loc='upper right')

    if 'DIR' in tau_inst.rev_info:
        # for chords, check direction order
        d0 = rho_km[1] - rho_km[0]
        d1 = rho_km[-1] - rho_km[-2]

        if d0 < 0:
            # ingress first, ind0 will be first index of occ within ring bounds
            ind0 = np.argwhere(rho_km < 150000.)[0][0]
            ind1 = np.argwhere(rho_km>150000.)[-1][0]
        else:
            ind0_1 = ind1
            ind1_1 = ind0
            ind0 = ind0_1
            ind1 = ind1_1

        ax1.scatter(rho_km[ind0]/1000., t_oet_1kspm[ind0], s=50,
                edgecolor='b', facecolor='none')
        ax1.scatter(rho_km[ind0]/1000., t_set_1kspm[ind0], s=50,
                edgecolor='r', facecolor='none')
    
        ax3.scatter(rho_km[ind0]/1000., phi_rl_deg[ind0], s=50,
                edgecolor='b', facecolor='none')
        ax3.scatter(rho_km[ind0]/1000., phi_ora_deg[ind0], s=50,
                edgecolor='r', facecolor='none')

        ax3.scatter(rho_km[ind0]/1000., phi_rl_deg[ind0], s=50,
                edgecolor='b', facecolor='none')
        ax3.scatter(rho_km[ind0]/1000., phi_rl_deg[ind0], s=50,
                edgecolor='b', facecolor='none')


        ax4.scatter(rho_km[ind0]/1000., phi_rl_dot_kms[ind0], s=50,
                edgecolor='b', facecolor='none')
        ax4.scatter(rho_km[ind0]/1000., rho_dot_kms[ind0], s=50,
                edgecolor='r', facecolor='none')

        ax5.scatter(rho_km[ind0]/1000., F_km[ind0], s=50,
                edgecolor='r', facecolor='none')

    #radii = [74.490, 91.983, 117.516, 122.052, 133.424, 136.774, 140.461]
    radii = [74.490, 91.983, 117.516, 122.052, 136.774, 140.461]
    for row in range(nrow):
        for col in range(ncol):
            #axes[row, col].legend(loc = "upper right")
            axes[row, col].set_xlim(70., 150.)
            axes[row, col].set_xlabel('Ring Radius (1000 km)')
            axes[row, col].tick_params(axis='both', direction='in', top=True, bottom=True, left=True, right=True)
            for radius in radii:
                axes[row, col].axvline(x=radius, color='k', linestyle=':',
                                linewidth=0.5)
    
    plt.tight_layout()
    pdf.savefig()
    plt.close()


    return pdf

def plot_cal_overview(pdf, cal_inst, dlp_inst):
    t_oet_spm  = cal_inst.t_oet_spm_vals
    F_offset_fit= np.interp(dlp_inst.t_oet_spm_vals, t_oet_spm, cal_inst.f_offset_fit_vals)
    F_offset = np.interp(t_oet_spm, cal_inst.f_spm, cal_inst.f_offset)
    F_sky_hz = cal_inst.f_sky_hz_vals
    P_free = cal_inst.p_free_vals

    t_oet_spm_dlp = dlp_inst.t_oet_spm_vals

    P_free_dlp = np.interp(t_oet_spm_dlp, t_oet_spm, P_free)
    P_obs_dlp = savgol_filter(dlp_inst.p_norm_vals * P_free_dlp, 61, 3)

    floor_fsky = int(np.nanmax(F_sky_hz/1.e6))

    ytitle1 = '$f_{sky}$ - ' + str(floor_fsky) + ' (MHz)'
    ytitle2 = '$f_{offset}$ (Hz)'
    ytitle3 = '$P_{free}$ (2x10$^{9}$)'
    xtitle = 'Observed Event Time, $t_{ERT}$ (1000 s)'
    xlim = [min(dlp_inst.t_oet_spm_vals)/1000., max(dlp_inst.t_oet_spm_vals)/1000.]

    nrow = 3
    ncol = 1
    fig3, axes3 = plt.subplots(nrow, ncol, figsize=(8.5, 11), sharex=True)
    fig3.subplots_adjust(hspace=0)
    axes3[0].plot(t_oet_spm/1000., (F_sky_hz-(floor_fsky*1.e6))*1.e-6, color='b')
    axes3[0].set_ylabel(ytitle1)
    axes3[0].set_xlim(xlim)
    axes3[0].tick_params(axis='both', direction='in', top=True, bottom=True,
            left=True, right=True)
    
    axes3[1].plot(t_oet_spm/1000., F_offset, color='b', linewidth=0.6)
    axes3[1].plot(dlp_inst.t_oet_spm_vals/1000., F_offset_fit, color='r',
            linewidth=0.6)
    axes3[1].set_ylabel(ytitle2)
    axes3[1].set_xlim(xlim)
    axes3[1].set_ylim([np.nanmin(F_offset_fit)-0.1, 
                        np.nanmax(F_offset_fit)+0.1])
    axes3[1].yaxis.tick_right()
    axes3[1].yaxis.set_label_position("right")
    axes3[1].tick_params(axis='both', direction='in', top=True, bottom=True,
            left=True, right=True)

    axes3[2].plot(t_oet_spm_dlp/1000., P_obs_dlp*2.e-9, color='b', linewidth=0.6)
    axes3[2].plot(t_oet_spm_dlp/1000., P_free_dlp*2.e-9, color='r', linewidth=0.6)
    axes3[2].tick_params(axis='both', direction='in', top=True, bottom=True,
            left=True, right=True)
    axes3[2].set_ylabel(ytitle3)
    axes3[2].set_xlabel(xtitle)
    axes3[2].set_xlim(xlim)
    axes3[2].set_ylim([0.,max(P_free_dlp*2.e-9)+0.2])
    pdf.savefig()
    plt.close()
    return pdf
def plot_tau_overview(pdf, geo_inst, tau_inst):
    tau = tau_inst.tau_vals
    tau_thresh = tau_inst.tau_threshold_vals
    rho_tau = tau_inst.rho_km_vals
    elev_deg = geo_inst.elev_deg_vals
    rho_geo = geo_inst.rho_km_vals


    tau = tau_inst.tau_vals
    tau1 = savgol_filter(tau, 123, 3)
    rho_tau = tau_inst.rho_km_vals
    elev_deg = np.interp(tau_inst.t_oet_spm_vals, geo_inst.t_oet_spm_vals,
            geo_inst.elev_deg_vals)
    rho_geo = geo_inst.rho_km_vals



    dsn = geo_inst.rev_info['dsn']

    #ylim1= [5.5, -0.19]
    max_thresh = np.nanmax(tau_thresh)
    if max(tau1) > 5.5:
        ylim1 = [5.5, -0.2]
    else:
        ylim1 = [max(tau1), -0.2]
    xlim1 = [70., 145.]

    ylim2 = [0., 90.]

    fig, ax1 = plt.subplots(figsize=(11,7))
    plt.title('Cassini RSS: Rev' + tau_inst.rev_info['rev_num'] + '-'
            + tau_inst.rev_info['prof_dir'][1:-1] + ', '
            + tau_inst.rev_info['band'][1:-1]
            + '-Band, ' + tau_inst.rev_info['dsn'])
    ax2 = ax1.twinx()
    ax1.tick_params(axis='x', direction='in', top=True, bottom=True)
    ax1.tick_params(axis='y', colors='blue', direction='in', right=True, left=True)
    ax2.tick_params(axis='y', colors='magenta', direction='in')

    ax1.yaxis.label.set_color('blue')
    ax2.yaxis.label.set_color('magenta')

    ax1.plot(rho_tau/1000., tau, color='blue', linewidth=0.3)
    ax1.plot(rho_tau/1000., tau_thresh, linestyle='--', color='red', linewidth=0.8)
    ax1.set_ylabel('Normal Optical Depth')
    ax1.set_xlabel('Ring Radius, $\\rho$ (1000 km)')
    ax1.spines['left'].set_color('blue')
    ax1.set_xlim(xlim1)
    ax1.set_ylim(ylim1)
    ax1.grid(b=True)

    ax2.plot(rho_tau/1000., elev_deg, color='magenta', linewidth=0.8)
    ax2.set_xlim(xlim1)
    ax2.set_ylabel(dsn + ' Elevation Angle (deg)')
    ax2.spines['right'].set_color('magenta')
    ax2.spines['left'].set_color('blue')
    ax2.set_ylim(ylim2)

    # Cring
    yt = -0.07
    ax1.text(82, yt, 'C')
    ax1.text(94, yt, 'B1')
    ax1.text(101, yt, 'B2')
    ax1.text(106, yt, 'B3')
    ax1.text(112, yt, 'B4')
    ax1.text(118, yt, 'CD')
    ax1.text(127, yt, 'A')
    pdf.savefig()

    plt.close()
    return pdf

def plot_tau(pdf, tau_inst):
    # Plot page 5-9 -- optical depth as a funct of ring radius
    rho_km = tau_inst.rho_km_vals
    tau = tau_inst.tau_vals
    tau_thresh = tau_inst.tau_threshold_vals
    res_km = str(round(tau_inst.input_res,3))
    band = str(tau_inst.rev_info['band'].split('"')[1])
    if band=='K':
        band='Ka'

    
    title = ('Cassini RSS: Reconstructed X-band Normal Optical Depth Profile'
            + ' (' + res_km + ' km Resolution)')
    ncol = 1
    nrow = 4

    
    rho_min = 74.
    rho_max = 155.

    ylim_max = round((max(tau_thresh) * 1.5)*2.0) / 2.0
    ylim = [ylim_max, -0.5]
    yticks_max = int(np.floor(ylim_max))
    yticks=range(0, yticks_max)

    
    for page in range(17):
        fig, ax = plt.subplots(nrow, ncol, figsize=(8.5,11))
        plt.setp(ax, yticks=yticks)#yticks=[0,1,2,3,4,5])
         #plt.locator_params(axis='y', nbins=6)
        for n in range(nrow):
            ax[n].axhline(y=0, color='c')
            ax[n].plot(rho_km/1000., tau_thresh, color='r', linestyle='--')
            ax[n].plot(rho_km/1000., tau, color='b', linewidth=0.8)
            ax[n].tick_params(axis='both', direction='in', top=True, bottom=True, left=True, right=True)
            r2 = rho_min + 1.
            ax[n].set_xlim(rho_min, r2)
            ax[n].set_ylim(ylim)
            ax[n].grid(True, which='both')
            ax[n].set_ylabel('$\\tau$')
            if n==0:
                ax[n].set_title(title)
            if n==(nrow-1):
                ax[n].set_xlabel('Ring Radius, $\\rho$ (1000 km)')
                
            
            rho_min = r2
        plt.tight_layout()
        pdf.savefig()
        plt.close()
    return pdf

def plot_phase(pdf, tau_inst):
    res_km = str(round(tau_inst.input_res,3))
    band = str(tau_inst.rev_info['band'].split('"')[1])
    if band=='K':
        band='Ka'


    # Plot page 10 -- phase shift as a funct of ring radius
    rho_km = tau_inst.rho_km_vals
    phi_deg = np.degrees(tau_inst.phase_vals)



    title = ('Cassini RSS: Reconstructed ' + band + '-band Phase Shift Profile '
                +'(' + res_km + ' km Resolution)')
    ytitle = '$\phi$ (deg.)'
    xtitle = 'Ring Radius, $\\rho$ (1000 km)'
    ncol = 1
    nrow = 4
    fig5, axes5 = plt.subplots(nrow, ncol, figsize=(8.5,11))

    axes5[0].set_title(title)
    axes5[0].plot(rho_km/1000., phi_deg, color='b', linewidth=lw1)
    axes5[0].set_xlim([72., 90.])
    axes5[0].set_ylabel(ytitle)
    axes5[0].tick_params(axis='both', direction='in', top=True, bottom=True, left=True, right=True)

    axes5[1].plot(rho_km/1000., phi_deg, color='b', linewidth=lw1)
    axes5[1].set_xlim([90., 108.])
    axes5[1].set_ylabel(ytitle)
    axes5[1].tick_params(axis='both', direction='in', top=True, bottom=True, left=True, right=True)

    axes5[2].plot(rho_km/1000., phi_deg, color='b', linewidth=lw1)
    axes5[2].set_xlim([108., 126.])
    axes5[2].set_ylabel(ytitle)
    axes5[2].tick_params(axis='both', direction='in', top=True, bottom=True, left=True, right=True)

    axes5[3].plot(rho_km/1000., phi_deg, color='b', linewidth=lw1)
    axes5[3].set_xlim([126., 144.])
    axes5[3].set_ylabel(ytitle)
    axes5[3].set_xlabel(xtitle)
    axes5[3].tick_params(axis='both', direction='in', top=True, bottom=True, left=True, right=True)
    pdf.savefig()
    plt.close()
    return pdf

def plot_summary_doc_v2(geo_inst, cal_inst, dlp_inst, tau_inst):
    pd1 = (tau_inst.rev_info['prof_dir'].split('"')[1])[0]
    if 'DIR' in tau_inst.rev_info.keys():
        pd1 = 'C' + pd1
    if 'PER' in tau_inst.rev_info.keys():
        pd1 = 'P' + pd1
    revstr = tau_inst.rev_info['rev_num'].zfill(3)
    outtitle, outdir = construct_filepath(tau_inst.rev_info, 'Summary')
    outfig = outdir[0] + 'Rev' + revstr + pd1 + '_' + outtitle[0] + '.pdf'
    with PdfPages(outfig) as pdf:
        pdf = plot_bullseye(pdf, dlp_inst)

        pdf = plot_occ_earth_view(pdf, geo_inst)
        pdf = plot_occ_pole_view(pdf, geo_inst)
        pdf = plot_geo_overview(pdf, geo_inst, tau_inst)
        pdf = plot_cal_overview(pdf, cal_inst, dlp_inst)
        pdf = plot_tau_overview(pdf, geo_inst, tau_inst)
        pdf = plot_tau(pdf, tau_inst)
        pdf = plot_phase(pdf, tau_inst)
    latex_summary_doc(outfig, tau_inst.input_res, outfig[:-4])
    print('\tSaving summary PDF to:\n\t\t' + outfig)


    return None





