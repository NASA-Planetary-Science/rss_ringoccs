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
 *  Purpose:                                                                  *
 *      Provides tools for modeling ideal diffraction data.                   *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       June 2, 2025                                                  *
 ******************************************************************************/

/*  Include guard to prevent including this file twice.                       */
#ifndef RSS_RINGOCCS_MODEL_H
#define RSS_RINGOCCS_MODEL_H

/*  TMPL_RESTRICT given here.                                                 */
#include <libtmpl/include/tmpl_config.h>

/*  size_t typedef provided here.                                             */
#include <stddef.h>

/*  Tau object is typedef'd here.                                             */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

enum rssringoccs_ModelType {
    rssringoccs_Model_LeftEdge,
    rssringoccs_Model_RightEdge,
    rssringoccs_Model_Ringlet,
    rssringoccs_Model_Gap,
    rssringoccs_Model_SquareWave,
    rssringoccs_Model_SineWave,
    rssringoccs_Model_DiracDelta,
    rssringoccs_Model_None
};

typedef struct rssringoccs_ModelParameters_Def {
    enum rssringoccs_ModelType model;
    double wavelength;
    double peak_opacity;
    union {
        struct {
            double center;
            double width;
        } square;
        struct {
            double center;
        } edge;
        struct {
            double left;
            double right;
            size_t number_of_waves;
        } wave;
    } geometry;
} rssringoccs_ModelParameters;

void
rssringoccs_Tau_Model_Left_Straightedge(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const rssringoccs_ModelParameters * TMPL_RESTRICT const parameters
);

void
rssringoccs_Tau_Model_Right_Straightedge(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const rssringoccs_ModelParameters * TMPL_RESTRICT const parameters
);

#endif
/*  End of include guard.                                                     */
