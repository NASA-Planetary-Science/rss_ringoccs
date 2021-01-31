/******************************************************************************
 *                                 LICENSE                                    *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify it   *
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
 *                          rss_ringoccs_linalg                               *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Define matrices and vectors and provide the prototypes for basic      *
 *      linear algebra functions and routines.                                *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       January 26, 2021                                              *
 ******************************************************************************/

/*  Include guard to prevent including this file twice.                       */
#ifndef __RSS_RINGOCCS_LINALG_H__
#define __RSS_RINGOCCS_LINALG_H__

#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

typedef struct rssringoccs_DoubleMatrix {
    double **data;
    unsigned long row_length;
    unsigned long column_length;
    rssringoccs_Bool error_occurred;
    char *error_message;
} rssringoccs_Matrix;

typedef struct rssringoccs_ComplexDoubleMatrix {
    rssringoccs_ComplexDouble **data;
    unsigned long row_length;
    unsigned long column_length;
    rssringoccs_Bool error_occurred;
    char *error_message;
} rssringoccs_ComplexDoubleMatrix;

#endif
/*  End of #ifndef __RSS_RINGOCCS_LINALG_H__.                                 */
