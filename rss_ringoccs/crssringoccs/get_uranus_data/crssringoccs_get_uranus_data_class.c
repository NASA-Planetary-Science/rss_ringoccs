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
 *      Defines the GetUranusData class for rss_ringoccs. This is used to     *
 *      extract geometry (Geo), calibration (Cal), diffraction-limited (DLP), *
 *      and reconstructed (Tau) data from CSV (or .TAB) files and create a    *
 *      Python object out of them. This way the data is more easily worked    *
 *      with as numpy arrays.                                                 *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       June 22, 2019                                                 *
 ******************************************************************************/
#include "../crssringoccs.h"

PyTypeObject crssringoccs_GetUranusData = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "GetUranusData",
    .tp_doc =
        "\r\tPurpose:\n"
        "\r\t\tExtracts data from the column separated TAB files and\n"
        "\r\t\tcreates numpy arrays out of them. All of the data is collected\n"
        "\r\t\tinto a single class, an instance of GetUranusData.\n"
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
        "\r\t\t\tdirectly with the GetUranusData class (by passing it\n"
        "\r\t\t\tto DiffractionCorrection) with reconstructions on the PDS.\n"
        "\r\t\t2.)\n"
        "\r\t\t\tThe geometry and calibration data is interpolated against\n"
        "\r\t\t\tthe DLP data. This assumes the DLP and CAL file represent\n"
        "\r\t\t\ta single occultation. For chord occultations, this means\n"
        "\r\t\t\tthe data must be separated into two portions, ingress and\n"
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
        "\r\t\t\tset use_deprecate = True.\n",
    .tp_basicsize = sizeof(crssringoccs_PyCSVObj),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = PyType_GenericNew,
    .tp_init = (initproc)crssringoccs_GetUranusData_Init,
    .tp_dealloc = (destructor)crssringoccs_PyCSVObj_Destroy,
    .tp_members = crssringoccs_PyCSVObj_Members,
    .tp_methods = crssringoccs_PyCSVObj_Methods
};
