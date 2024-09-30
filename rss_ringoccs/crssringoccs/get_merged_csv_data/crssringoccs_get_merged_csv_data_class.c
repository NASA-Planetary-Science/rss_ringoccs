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
 *                   crssringoccs_get_merged_csv_data_class                   *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Defines the GetMergedCSVData class for rss_ringoccs. This is used to  *
 *      extract data from the merged DLP files (DLPM.TAB) and convert them    *
 *      into a Python object (the GetMergedCSVData class).                    *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       June 22, 2019                                                 *
 ******************************************************************************/
#include "../crssringoccs.h"

PyTypeObject crssringoccs_GetMergedCSVData = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "GetMergedCSVData",
    .tp_doc =
        "\r\tPurpose:\n"
        "\r\t\tExtracts data from a column separated DLPM.TAB file and\n"
        "\r\t\tcreates a Python object consisting of each column as an\n"
        "\r\tattribute. All of the data is collected\n"
        "\r\t\tinto a single class, an instance of GetUranusData.\n"
        "\r\tArguments:\n"
        "\r\t\tdlpm (str):\n"
        "\r\t\t\tPath to the merged DLP file (DLPM.TAB).\n"
        "\r\tNotes:\n"
        "\r\t\t1.)\n"
        "\r\t\t\tDLPM CSV files should have all of the data already\n"
        "\r\t\t\tinterpolated against rho_km_vals. The variables should be\n"
        "\r\t\t\tfunctions of rho_km_vals as well.\n",
    .tp_basicsize = sizeof(crssringoccs_PyCSVObj),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = PyType_GenericNew,
    .tp_init = (initproc)crssringoccs_GetMergedCSVData_Init,
    .tp_dealloc = (destructor)crssringoccs_PyCSVObj_Destroy,
    .tp_members = crssringoccs_PyCSVObj_Members,
    .tp_methods = crssringoccs_PyCSVObj_Methods
};
