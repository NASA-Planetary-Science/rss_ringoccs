
"""
plot_rss_occ.py

Purpose: Plot an earth-view of a Saturn ring occultation.

Revisions:
    2018 Jul 24 - jfong - copied from plot_rss_occ_v3.py
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pds3_reader import PDS3Reader
import pdb
import time
import spiceypy as spice

radii = [60268., 60268., 54364.]
rings_km = [74490., 91983., 117516., 122052., 136774., 139826.]

def plot_rss_occ(pdf, B_deg_mean, phi_ora_deg_vals, rho_km_vals, 
        t_oet_spm_vals):

    # Initialize 3D subplot with no axes or gridlines
    az = 0.
    elev = B_deg_mean
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.set_axis_off()
    ax.grid(b=False)
    ax.set_xlim([-1.e5, 1.e5])
    ax.set_ylim([-1.e5, 1.e5])
    ax.set_zlim([-1.e5, 1.e5])

    # Plot saturn ellipsoid on both plots
    u = np.linspace(0, 2 * np.pi, 25)
    v = np.linspace(0, np.pi, 16)
    
    rx = radii[0] * np.outer(np.cos(u), np.sin(v))
    ry = radii[1] * np.outer(np.sin(u), np.sin(v))
    rz = radii[2] * np.outer(np.ones_like(u), np.cos(v))

    surf = ax.plot_surface(rx, ry, rz, linewidth=0.5, edgecolors='black',
                        color='white', shade=False, alpha=1)
    # Plot rings
    ntheta = 361
    theta = np.linspace(0., np.pi*2., ntheta)
    cosTheta = np.cos(theta)
    sinTheta = np.sin(theta)

    Brad = np.radians(B_deg_mean)
    Re = radii[1]
    Rp = radii[2]
    Rpsky = np.sqrt((Re * np.sin(Brad))**2 + (Rp*np.cos(Brad))**2)
    rot_ax = [0., 1., 0.]
    rot_mat = spice.axisar(rot_ax, (Brad))

    ax.view_init(elev, az)

    # Plot x (yellow), y (red), and z(green) axes
#    ax.plot3D([0, 100000.], [0, 0], [0, 0], 'y-', linewidth=2.3)
#    ax.plot3D([0, 0], [0, 100000.], [0, 0], 'r-', linewidth=2.3)
#    ax.plot3D([0, 0], [0, 0], [0, 100000.], 'g-', linewidth=2.3)

    
    # Plot rings 
    for nr in rings_km:
        xr = nr * cosTheta
        yr = nr * sinTheta
        zr = np.zeros(len(xr))
        p_ind = ind_blocked_by_saturn(xr, yr, zr, rot_mat, Rpsky, Re, 
                log_op='or')
        yvals = yr[p_ind] 
        xvals = -xr[p_ind] 
        zvals = zr[p_ind] 


        
        ax.plot3D(xvals, yvals, zvals, color='black', linewidth=0.5)

    # Plot occultation track
    occx = rho_km_vals * np.cos(np.radians(phi_ora_deg_vals))
    occy = rho_km_vals * np.sin(np.radians(phi_ora_deg_vals))
    occz = np.zeros(len(occx))


    # check if occ blocked by planet, if so use blue dashes
    p_ind_occ = ind_blocked_by_saturn(occx, occy, occz, rot_mat, Rpsky, Re)



    ax.plot3D(occx[p_ind_occ], occy[p_ind_occ], color='blue')
    ax.plot3D(occx[p_ind_occ==False], occy[p_ind_occ==False], color='blue',
            linestyle='--')

    # add 30-min ticks on occ trqck
    #   NOTE: this assumes 1sec intervals between points
    dpts = 60 * 30
    ax.scatter(occx[::dpts], occy[::dpts], color='blue', linewidth=1.5,
            alpha=1)

    # Plot marker in direction of earth
    rlen_km2 = 1.6e5

    ax.plot3D([68000., rlen_km2+20000.], [0, 0], [0, 0], 'k-', linewidth=2.3)
    pdf.savefig()
    return pdf

def ind_blocked_by_saturn(xvals, yvals, zvals, rot_mat, Rpsky, Re, log_op='and'):
    npts = len(xvals)
    xsky = np.zeros(npts)
    ysky = np.zeros(npts)
    zsky = np.zeros(npts)
    for n in range(npts):
        vec = [xvals[n], yvals[n], zvals[n]]
        vec_rot = np.dot(rot_mat, vec)
        xsky[n] = vec_rot[1]
        ysky[n] = vec_rot[2]
        zsky[n] = vec_rot[0]
    Rsky = np.sqrt(xsky**2 + (Re*ysky/Rpsky)**2)
    if log_op == 'and':
        p_ind = np.logical_and(Rsky>Re, zsky<0)
    else:
        p_ind = np.logical_or(Rsky>Re, zsky<0)

    return p_ind



def main():
    # Read geometry data file
    data_dir = '../../../data/Cassini_RSS_Ring_Profiles_2018_Archive/Rev7E/Rev007E_RSS_2005_123_X43_E/'
    geo_tab = data_dir + 'RSS_2005_123_X43_E_GEO.TAB'
    geo_lbl = data_dir + 'RSS_2005_123_X43_E_GEO.LBL'

    geo_lbl_inst = PDS3Reader(geo_lbl)
    B_deg_mean = np.mean(geo_lbl_inst.series.B_deg_vals)
#    B_deg_mean=90.
    phi_ora_deg_vals = geo_lbl_inst.series.phi_ora_deg_vals
    rho_km_vals = geo_lbl_inst.series.rho_km_vals
    t_oet_spm_vals = geo_lbl_inst.series.t_oet_spm_vals
    plot_rss_occ(B_deg_mean, phi_ora_deg_vals, rho_km_vals,
            t_oet_spm_vals, 'test_rss_occ', SAVE=0)
    pdb.set_trace()


if __name__ == "__main__":
    main()


