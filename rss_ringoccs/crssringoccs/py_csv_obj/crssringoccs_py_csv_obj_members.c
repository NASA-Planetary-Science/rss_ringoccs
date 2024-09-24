/******************************************************************************
 *                                  LICENSE                                   *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  rss_ringoccs is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************/
#include "../crssringoccs.h"

PyMemberDef crssringoccs_PyCSVObj_Members[] = {

    /*  Ring radius, primary independent variable for the all of the data.    */
    {
        "rho_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, rho_km_vals),
        0,
        "Ring radius, in kilometers. This is the radial distance from\n"
        "the core of the planet to the point of interest."
    },

    /*  Complex phase angle for the diffracted transmittance.                 */
    {
        "phase_deg_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, phase_deg_vals),
        0,
        "Raw diffracted phase, in degrees. The diffracted transmittance\n"
        "is a point in the complex plane. This variable represents the angle\n"
        "this point makes with the positive real axis. It varies\n"
        "as a function of the ring radius, rho_km_vals."
    },

    /*  Opening angle for planetary (i.e. Saturn) rings with respect to Earth.*/
    {
        "B_deg_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, B_deg_vals),
        0,
        "Ring inclination angle, in degrees. Also called the opening angle.\n"
        "This is the angle the ring plane makes with the line going from\n"
        "Earth to the spacecraft (i.e. Cassini)."
    },

    /*  Spacecraft-ring distance.                                             */
    {
        "D_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, D_km_vals),
        0,
        "Spacecraft to ring-intercept point distance, in kilometers. This\n"
        "is the distance between the point in the plane the data point\n"
        "represents, and position of the observer at the time of\n"
        "observation. Note that it can take up to a second for light\n"
        "to travel from the observer (i.e. Cassini) to the ring plane.\n"
        "Because of this, the location of the observation at the time of\n"
        "observation and the location of the observer at the event time\n"
        "(when the light reaches the rings) can differ significantly."
    },

    /*  x-component of the spacecraft.                                        */
    {
        "rx_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, rx_km_vals),
        0,
        "x coordinate of the spacecraft in the planetocentric frame,\n"
        "in kilometers. The spacecraft lies at the point\n"
        "(rx_km_vals, ry_km_val, rz_km_vals) in space. The planetocentric\n"
        "frame takes the core of the planet as the center and the rotation\n"
        "axis as the z-axis. The x and y axis are chosen using a right-hand\n"
        "rule convention. See MTR86 for details."
    },

    /*  y-component of the spacecraft.                                        */
    {
        "ry_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, ry_km_vals),
        0,
        "y coordinate of the spacecraft in the planetocentric frame,\n"
        "in kilometers. The spacecraft lies at the point\n"
        "(rx_km_vals, ry_km_val, rz_km_vals) in space. The planetocentric\n"
        "frame takes the core of the planet as the center and the rotation\n"
        "axis as the z-axis. The x and y axis are chosen using a right-hand\n"
        "rule convention. See MTR86 for details."
    },

    /*  y-component of the spacecraft.                                        */
    {
        "rz_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, rz_km_vals),
        0,
        "y coordinate of the spacecraft in the planetocentric frame,\n"
        "in kilometers. The spacecraft lies at the point\n"
        "(rx_km_vals, ry_km_val, rz_km_vals) in space. The planetocentric\n"
        "frame takes the core of the planet as the center and the rotation\n"
        "axis as the z-axis. The x and y axis are chosen using a right-hand\n"
        "rule convention. See MTR86 for details."
    },

    /*  Frequency of the incoming radio signal.                               */
    {
        "f_sky_hz_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, f_sky_hz_vals),
        0,
        "Frequency of the input signal, in Hertz.\n"
        "This is computed as the difference of the predicted\n"
        "frequency and the residual frequency."
    },

    /*  Diffracted power.                                                     */
    {
        "p_norm_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, p_norm_vals),
        0,
        "Raw power data, unitless. This has been normalized so that\n"
        "the free-space region has power 1. Since the wavelength of\n"
        "the incoming light is so large, the planetary rings cause\n"
        "diffraction to occur. The measured signal is thus the diffracted\n"
        "power. p_norm_vals contains this data, and not the diffraction\n"
        "corrected profile."
    },

    /*  Diffracted optical depth.                                             */
    {
        "raw_tau_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, raw_tau_vals),
        0,
        "Raw optical depth, unitless. This is scaled logarithm of the\n"
        "diffracted power p_norm_vals, where power is scaled by\n"
        "a geometric factor involving the ring opening angle."
    },

    /*  Polar angle for the point of interest.                                */
    {
        "phi_deg_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, phi_deg_vals),
        0,
        "Ring azimuth angle, in degrees. This is the polar angle of the\n"
        "point of interest. The ordered pair (rho_km_vals, phi_deg_vals)\n"
        "uniquely specifies the location of the ring-intercept point\n"
        "in the ring plane."
    },

    /*  Longitudinal angle.                                                   */
    {
        "phi_rl_deg_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, phi_rl_deg_vals),
        0,
        "Ring longitude angle, in degrees."
    },

    {
        "raw_tau_threshold_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, raw_tau_threshold_vals),
        0,
        "Raw threshold optical depth. Unitless."
    },

    {
        "rev_info",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, rev_info),
        0,
        "Information about the occultation. This is set to 'None' since\n"
        "the ExtractCSVData class is not provided an RSR instance to fetch\n"
        "the rev info from. It is kept for the sake of compatibility with\n"
        "other classes (namely DiffractionCorrection)."
    },

    {
        "rho_corr_pole_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, rho_corr_pole_km_vals),
        0,
        "Ring radius with pole correction, in kilometers."
    },

    {
        "rho_corr_timing_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, rho_corr_timing_km_vals),
        0,
        "Ring radius with timing correction, in kilometers."
    },

    {
        "rho_dot_kms_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, rho_dot_kms_vals),
        0,
        "Time derivative of the ring radius, in kilometers per second.\n"
        "This is the radial time derivative, the amount of change in\n"
        "rho over a small amount of time. It is kept as a function of rho\n"
        "and not as a function of time."
    },

    {
        "t_oet_spm_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, t_oet_spm_vals),
        0,
        "Observed event time in seconds past midnight. This is the\n"
        "time the signal was received on Earth."
    },

    {
        "t_ret_spm_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, t_ret_spm_vals),
        0,
        "Ring event time in seconds past midnight. This is the moment of\n"
        "time when the signal crossed the ring plane."
    },

    {
        "t_set_spm_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, t_set_spm_vals),
        0,
        "Spacecraft event time in seconds past midnight. This is when the\n"
        "signal left the spacecraft. For many data sets with Cassini and\n"
        "Saturn there is about a 1 light second distance between the\n"
        "spacecraft and the ring intercept point, meaning\n"
        "t_ret_spm_vals - t_set_spm_vals ~ 1 second."
    },

    {
        "tau_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, tau_vals),
        0,
        "Corrected optical depth, unitless. Set to 'None' if a tau file\n"
        "is not provided. This is the scaled logarithmic power, but with\n"
        "diffraction correction applied."
    },

    {
        "tau_power_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, tau_power_vals),
        0,
        "Normalized optical power, unitless and corrected for diffraction\n"
        "effects. This is the output of the Fresnel inversion process on\n"
        "the diffracted data."
    },

    {
        "tau_phase_deg_vals",
        T_OBJECT_EX, offsetof(crssringoccs_PyCSVObj, tau_phase_deg_vals),
        0,
        "Diffraction corrected phase, in degrees. The diffraction corrected\n"
        "complex transmittance is a point in the complex plane.\n"
        "This is the angle, in degrees, the point makes with the\n"
        "real axis."
    },

    {
        "history",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, history),
        0,
        "History of the object. This includes the following data:\n"
        "rss_ringoccs Version.\n"
        "libtmpl Version.\n"
        "Python Version.\n"
        "Host Name.\n"
        "User Name.\n"
        "Run Date.\n"
        "Operating System.\n"
        "Positional Arguments (geo, cal, dlp).\n"
        "Keyword Arguments (tau, use_deprecated).\n"
    },

    /*  Terminator for the class.                                             */
    {
        NULL
    }
};
