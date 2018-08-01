"""

create_summary_doc.py

Purpose: Create the plots in the second version of Essam's EASYDATA summary PDF.

Revisions:
    2018 Jul 23 - jfong - copied from create_summary_doc_v2.py
"""


import numpy as np
import matplotlib
from pds3_reader import PDS3Reader
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pdb
import time
from plot_rss_occ import plot_rss_occ

sat_radius =  60268.
rings_km = [74490., 91983., 117516., 122052., 136774., 139826.]
radii = [74.490, 91.983, 117.516, 122.052, 133.424, 136.774, 140.461]

def plot_summary_doc(geo_inst, cal_inst, dlp_inst, tau_inst, outfig):
    # Extract relevant informatio from instances
    rho_km = geo_inst.rho_km_vals
    phi_rl_rad = np.radians(geo_inst.phi_rl_deg_vals)
    phi_ora_rad = np.radians(geo_inst.phi_ora_deg_vals)
    B_mean = np.mean(geo_inst.B_deg_vals)
    band = geo_inst.FREQUENCY_BAND #NOTE: hardcoded in because not avail in instances read from lbl
    t_oet_spm_geo = geo_inst.t_oet_spm_vals
    t_oet_hrs = geo_inst.t_oet_spm_vals / 60. / 60.
    t_set_hrs = geo_inst.t_set_spm_vals / 60. / 60.

    D_km = geo_inst.D_km_vals
    phi_rl_deg = geo_inst.phi_rl_deg_vals
    phi_ora_deg = geo_inst.phi_ora_deg_vals
    
    rho_dot_kms = geo_inst.rho_dot_kms_vals
    phi_rl_dot_kms = geo_inst.phi_rl_dot_kms_vals
    F_km = geo_inst.F_km_vals
    R_imp_km = geo_inst.R_imp_km_vals

    t_oet_spm = cal_inst.t_oet_spm_vals
    F_sky_hz = cal_inst.f_sky_hz_vals
    F_resid_hz = cal_inst.f_sky_resid_fit_vals
    P_free = cal_inst.p_free_vals
    
    t_oet_spm_dlp = dlp_inst.t_oet_spm_vals
    try:
        P_norm_dlp = dlp_inst.p_norm_vals
    except AttributeError:
        mu = np.sin(abs(np.radians(dlp_inst.B_deg_vals)))
        P_norm_dlp = np.exp(-dlp_inst.tau_vals/mu)


    P_free_interp = np.interp(t_oet_spm_dlp, t_oet_spm, P_free)
    P_obs_dlp = P_norm_dlp * P_free_interp

    with PdfPages(outfig) as pdf:
        # Plot cover page -- bird's eye view of ring plane with occ tracks
#        rho_km = geo_inst.rho_km_vals
#        phi_rl_rad = np.radians(geo_inst.phi_rl_deg_vals)
#        phi_ora_rad = np.radians(geo_inst.phi_ora_deg_vals)
#        B_mean = np.mean(geo_inst.B_deg_vals)
        
        occ_ed_x = rho_km * np.cos(phi_ora_rad)
        occ_ed_y = rho_km * np.sin(phi_ora_rad)
        
        occ_an_x = rho_km * np.cos(phi_rl_rad)
        occ_an_y = rho_km * np.sin(phi_rl_rad)
        
        theta = np.linspace(0., 2.*np.pi, 361)
        plt.figure()

        xx = np.cos(theta)
        yy = np.sin(theta)

        sat_xx = sat_radius/1000. * xx
        sat_yy = sat_radius/1000. * yy


        plt.plot(sat_xx, sat_yy, 'k-', linewidth=1.5)

        plt.plot(occ_an_x/1000., occ_an_y/1000., 'b-', linewidth=1.5)
        plt.plot(occ_ed_x/1000., occ_ed_y/1000., 'r-', linewidth=1.5)
        plt.axis('equal')
        plt.xlim([-150., 150.])
        plt.ylim([-150., 150.])
        plt.xlabel('$x$ (1000 km)')
        plt.ylabel('$y$ (1000 km)')
        plt.title('Mean Observed Ring Elevation $B$ = '+str(B_mean)+'$^\circ$')
        plt.plot(0,0,'+')
        plt.text(5,-3, 'Saturn')
        plt.text(80, -3, 'C')
        plt.text(101, -3, 'B')
        plt.text(125, -3, 'A')

        for ring in rings_km:
            plt.plot(ring/1000.*xx, ring/1000.*yy, 'm-', linewidth=1.)

        pdf.savefig()
        plt.close()
        

        # Plot page 2  -- Earth-view and north-pole view of Saturn + occ tracks
        pdf = plot_rss_occ(pdf, B_mean, phi_ora_deg, rho_km, t_oet_spm_geo)
        plt.close()
        pdf = plot_rss_occ(pdf, 90., phi_ora_deg, rho_km, t_oet_spm_geo)
        plt.close()
        # Plot page 3 -- geometry params as a function of ring radius
        nrow = 3
        ncol = 2
#        fig, axes = plt.subplots(nrow, ncol, figsize=(7.5, 9))
        fig, axes = plt.subplots(nrow, ncol, figsize=(8.5, 11))
        ax1 = axes[0,0]
        ax2 = axes[0,1]
        ax3 = axes[1,0]
        ax4 = axes[1,1]
        ax5 = axes[2,0]
        ax6 = axes[2,1]

        # gridspec_kw = {'width_ratios':[3, 1]}
        
        
        ax1.set_title('Time (HPM 2005-123)', fontweight='bold', fontsize=10)
        ax1.plot(rho_km/1000., t_oet_hrs, 'b', label = 'ERT')
        ax1.plot(rho_km/1000., t_set_hrs, 'r', label='SCET')
        
        ax2.set_title('S/C Distance D/RS', fontweight='bold', fontsize=10)
        ax2.plot(rho_km/1000., D_km/sat_radius, 'k')
        
        ax3.set_title('Longitude, Azimuth (deg)',fontweight='bold', fontsize=10)
        ax3.plot(rho_km/1000., phi_rl_deg, 'b', label = '$\phi_E$')
        ax3.plot(rho_km/1000., phi_ora_deg, 'r', label = '$\phi_{J2K}$')
        
        
        ax4.set_title('Occ Track Velocity (km/s)', fontweight='bold', 
                 fontsize=10)
        ax4.plot(rho_km/1000., rho_dot_kms, 'r', label = '$V_{rad}$')
        ax4.plot(rho_km/1000., phi_rl_dot_kms, 'b', label = '$V_{az}$')
        
        ax5.set_title('Fresnel Scale (km)', fontweight='bold', fontsize=10)
        ax5.plot(rho_km/1000., F_km, 'r', label = str(band))
        
#        axes[2, 1].plot(rho_km/1000., R_imp_km, 'r', label = 'filler')
        ax6.set_title('Optical Depth Threshold', fontweight='bold',
                fontsize=10)
        
        radii = [74.490, 91.983, 117.516, 122.052, 133.424, 136.774, 140.461]
        for row in range(nrow):
            for col in range(ncol):
                axes[row, col].legend(loc = "upper right")
                axes[row, col].set_xlim(70., 150.)
                axes[row, col].set_xlabel('Ring Radius (1000 km)')
                for radius in radii:
                    axes[row, col].axvline(x=radius, color='k', linestyle=':',
                                    linewidth=0.5)
        
        plt.tight_layout()
        pdf.savefig()
        plt.close()

        # Plot page 4 -- calibration data
        ytitle1 = '$f_{sky}$ - 8426 (MHz)'
        ytitle2 = '$f_{resid}$ (Hz)'
        ytitle3 = '$P_{free}$'
        xtitle = 'Observed Event Time, $t_{ERT}$ (1000 s)'



        nrow = 3
        ncol = 1
        fig3, axes3 = plt.subplots(nrow, ncol, figsize=(8.5, 11), sharex=True)
        fig3.subplots_adjust(hspace=0)
        axes3[0].plot(t_oet_spm/1000., (F_sky_hz-8426.e6)*1.e-6, color='b')
        axes3[0].set_ylabel(ytitle1)
        
        axes3[1].plot(t_oet_spm/1000., F_resid_hz, color='r')
        axes3[1].set_ylabel(ytitle2)

        axes3[2].plot(t_oet_spm/1000., P_free, color='r')
        axes3[2].plot(t_oet_spm_dlp/1000., P_obs_dlp, color='g')
        axes3[2].set_ylabel(ytitle3)
        axes3[2].set_xlabel(xtitle)
        pdf.savefig()
        plt.close()

        # Plot page 5-9 -- optical depth as a funct of ring radius
        rho_km = tau_inst.rho_km_vals
        tau = tau_inst.tau_vals
        tau_thresh = tau_inst.tau_threshold_vals

        
        title = 'Cassini RSS: Reconstructed X-band Normal Optical Depth Profile (1 km Resolution)'
        ncol = 1
        nrow = 4

        
        rho_min = 74.
        rho_max = 155.

        ylim = [5.5, -0.5]

        
        for page in range(17):
            fig, ax = plt.subplots(nrow, ncol, figsize=(8.5,11))
            plt.setp(ax, yticks=[0,1,2,3,4,5])
#            plt.locator_params(axis='y', nbins=6)
            for n in range(nrow):
                ax[n].axhline(y=0, color='c')
#                ax[n].plot(rho_km/1000., tau_thresh, color='r')
                ax[n].plot(rho_km/1000., tau, color='b')
                r2 = rho_min + 1.
                ax[n].set_xlim(rho_min, r2)
                ax[n].set_ylim(ylim)
                ax[n].grid(True)
                ax[n].set_ylabel('$\\tau$')
                if n==0:
                    ax[n].set_title(title)
                if n==(nrow-1):
                    ax[n].set_xlabel('Ring Radius, $\\rho$ (1000 km)')
                    
                
                rho_min = r2
            plt.tight_layout()
            pdf.savefig()
            plt.close()
        # Plot page 10 -- phase shift as a funct of ring radius
        rho_km = tau_inst.rho_km_vals
        try:
            phi_deg = np.degrees(tau_inst.phase_vals)
        except AttributeError:
            phi_deg = tau_inst.phase_deg_vals


        title = 'Cassini RSS: Reconstructed X-band Phase Shift Profile (1 km Resolution)'
        ytitle = '$\phi$ (deg.)'
        xtitle = 'Ring Radius, $\\rho$ (1000 km)'
        ncol = 1
        nrow = 4
        fig5, axes5 = plt.subplots(nrow, ncol, figsize=(8.5,11))

        axes5[0].set_title(title)
        axes5[0].plot(rho_km/1000., phi_deg, color='b')
        axes5[0].set_xlim([72., 90.])
        axes5[0].set_ylabel(ytitle)

        axes5[1].plot(rho_km/1000., phi_deg, color='b')
        axes5[1].set_xlim([90., 108.])
        axes5[1].set_ylabel(ytitle)

        axes5[2].plot(rho_km/1000., phi_deg, color='b')
        axes5[2].set_xlim([108., 126.])
        axes5[2].set_ylabel(ytitle)

        axes5[3].plot(rho_km/1000., phi_deg, color='b')
        axes5[3].set_xlim([126., 144.])
        axes5[3].set_ylabel(ytitle)
        axes5[3].set_xlabel(xtitle)
        pdf.savefig()
        plt.close()
    return None


    




def main():
    # ***** Begin user input *****
    data_dir = '../../../../TC2017/data/Cassini_RSS_Ring_Profiles_2018_Archive/Rev7E/Rev007E_RSS_2005_123_X43_E/'
    
    geolblfile = data_dir + 'RSS_2005_123_X43_E_GEO.LBL'
    callblfile = data_dir + 'RSS_2005_123_X43_E_CAL.LBL'
    dlplblfile = data_dir + 'RSS_2005_123_X43_E_DLP_500M.LBL'
    taulblfile = data_dir + 'RSS_2005_123_X43_E_TAU_01KM.LBL'

    # Read label and data files
    geolbl_inst = PDS3Reader(geolblfile)
    callbl_inst = PDS3Reader(callblfile)
    dlplbl_inst = PDS3Reader(dlplblfile)
    taulbl_inst = PDS3Reader(taulblfile)
    current_time = time.strftime("_%Y%m%d-%H%M%S")
    outfig = 'plot_summary_doc'+current_time+'.pdf'
    plot_summary_doc(geolbl_inst.series, callbl_inst.series, dlplbl_inst.series, taulbl_inst.series, outfig)
    pdb.set_trace()
    
# ***** End user input *****
if __name__ == '__main__':
    main()




# Plot view from earth

