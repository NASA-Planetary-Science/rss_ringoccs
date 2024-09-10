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
 ******************************************************************************
 *                    crssringoccs_extract_csv_data_class                     *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Defines the ExtractCSVData class for rss_ringoccs. This is used to    *
 *      extract geometry (Geo), calibration (Cal), diffraction-limited (DLP), *
 *      and reconstructed (Tau) data from CSV (or .TAB) files and create a    *
 *      Python object out of them. This way the data is more easily worked    *
 *      with as numpy arrays.                                                 *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       June 22, 2019                                                 *
 ******************************************************************************/

/*  Data typoes and function prototypes provided here.                        */
#include "crssringoccs_extract_csv_data.h"


static PyMemberDef crssringoccs_ExtractCSVData_Members[] = {
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
        "is a point in the complex plane. This represents the angle\n"
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
        "represents, and position of the observation at the time of\n"
        "observation. Note that it can take up to a second for light\n"
        "to travel from the observer (i.e. Cassini) to the ring plane.\n"
        "Because of this, the location of the observation at the time of\n"
        "observation and the location of the observer at the event time\n"
        "(when the light reaches the rings) can differ significantly."
    },
    {
        "rx_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, rx_km_vals),
        0,
        "x coordinate of the spacecraft in planetocentric frame."
    },
    {
        "ry_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, ry_km_vals),
        0,
        "y coordinate of the spacecraft in planetocentric frame."
    },
    {
        "rz_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, rz_km_vals),
        0,
        "z coordinate of the spacecraft in planetocentric frame."
    },
    {
        "f_sky_hz_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, f_sky_hz_vals),
        0,
        "Frequency of the input signal"
    },
    {
        "p_norm_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, p_norm_vals),
        0,
        "Raw power data"
    },
    {
        "raw_tau_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, raw_tau_vals),
        0,
        "Raw optical depth"
    },
    {
        "phi_deg_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, phi_deg_vals),
        0,
        "Ring azimuth angle"
    },
    {
        "phi_rl_deg_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, phi_rl_deg_vals),
        0,
        "Ring longitude angle"
    },
    {
        "raw_tau_threshold_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, raw_tau_threshold_vals),
        0,
        "Raw threshold optical depth"
    },
    {
        "rev_info",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, rev_info),
        0,
        "Information about the occultation"
    },
    {
        "rho_corr_pole_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, rho_corr_pole_km_vals),
        0,
        "Ring radius with pole correction."
    },
    {
        "rho_corr_timing_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, rho_corr_timing_km_vals),
        0,
        "Ring radius with timing correction."
    },
    {
        "rho_dot_kms_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, rho_dot_kms_vals),
        0,
        "Time derivative of the ring radius."
    },
    {
        "t_oet_spm_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, t_oet_spm_vals),
        0,
        "Observed event time in seconds past midnight"
    },
    {
        "t_ret_spm_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, t_ret_spm_vals),
        0,
        "Ring event time in seconds past midnight"
    },
    {
        "t_set_spm_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, t_set_spm_vals),
        0,
        "Spacecraft event time in seconds past midnight"
    },
    {
        "tau_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, tau_vals),
        0,
        "Optical depth"
    },
    {
        "tau_power",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, tau_power),
        0,
        "Optical power"
    },
    {
        "tau_phase",
        T_OBJECT_EX, offsetof(crssringoccs_PyCSVObj, tau_phase),
        0,
        "Optical phase"
    },
    {
        "history",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, history),
        0,
        "History of the tau instance"
    },
    {
        "input_vars",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, input_vars),
        0,
        "Dictionary of input arguments used to create this instance."
    },
    {
        "input_kwds",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyCSVObj, input_kwds),
        0,
        "Dictionary of input keywords used to create this instance."
    },
    {
        NULL
    }
};

static PyMethodDef crssringoccs_ExtractCSVData_Methods[] =
{
    {
        NULL
    }
};

PyTypeObject ExtractCSVDataType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "ExtractCSVData",
    .tp_doc =
        "\r\tPurpose:\n"
        "\r\t\tExtracts data from the column separated TAB files and\n"
        "\r\t\tcreates numpy arrays out of them. All of the data is collected\n"
        "\r\t\tinto a single class, an instance of ExtractCSVData.\n"
        "\r\tArguments:\n"
        "\r\t\tgeo (str):\n"
        "\r\t\t\tPath to the geometry file (GEO.TAB).\n"
        "\r\t\tcal (str):\n"
        "\r\t\t\tPath to the calibration file (CAL.TAB).\n"
        "\r\t\tdlp (str):\n"
        "\r\t\t\tPath to the DLP file (DLP.TAB).\n"
        "\r\tKeywords:\n"
        "\r\t\ttau (str):\n"
        "\r\t\t\tPath to the reconstruction file (TAU.TAB).\n"
        "\r\t\tuse_deprecate (bool):\n"
        "\r\t\t\tBoolean for using the older file format from the PDS.\n"
        "\r\t\t\tTAB files that were generated before 2018 require this\n"
        "\r\t\t\tkeyword to be set to True.\n"
        "\r\tNotes:\n"
        "\r\t\t1.)\n"
        "\r\t\t\tThe tau file is not required. If provided, its values\n"
        "\r\t\t\tare interpolated against the values in the DLP file.\n"
        "\r\t\t\tThis is useful for comparing reconstructions made\n"
        "\r\t\t\tdirectly with the ExtractCSVData class (by passing it\n"
        "\r\t\t\tto DiffractionCorrection) with reconstructions on the PDS.\n"
        "\r\t\t2.)\n"
        "\r\t\t\tThe geometry and calibration data is interpolated against\n"
        "\r\t\t\tthe DLP data. This assumes the DLP and CAL file represent\n"
        "\r\t\t\ta single occultation. For chord occultations, this means\n"
        "\r\t\t\tthe data must be separated into two portion, ingress and\n"
        "\r\t\t\tegress. For GEO files the data does not need to be\n"
        "\r\t\t\tand the DLP file will be used to determine the occultation\n"
        "\r\t\t\ttype. The corresponding data will be extracted from GEO.\n"
        "\r\t\t3.)\n"
        "\r\t\t\tAll variables are made as functions of rho_km_vals.\n"
        "\r\t\t\tThat is, the ring radius is the independent variable among\n"
        "\r\t\t\tall of the data.\n"
        "\r\t\t4.)\n"
        "\r\t\t\tBy default, it is assumed the new formats are being used.\n"
        "\r\t\t\tThat is, use_deprecate is set to False as a default\n"
        "\r\t\t\tparameter. If older data sets are used (pre-2018), you must\n"
        "\r\t\t\tset use_deprecate = False.\n",
    .tp_basicsize = sizeof(crssringoccs_PyCSVObj),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = PyType_GenericNew,
    .tp_init = (initproc)crssringoccs_ExtractCSVData_Init,
    .tp_dealloc = (destructor)crssringoccs_ExtractCSVData_Destroy,
    .tp_members = crssringoccs_ExtractCSVData_Members,
    .tp_methods = crssringoccs_ExtractCSVData_Methods
};
