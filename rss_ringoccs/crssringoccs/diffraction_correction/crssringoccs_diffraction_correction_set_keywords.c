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

void
crssringoccs_DiffractionCorrection_Set_Keywords(
    crssringoccs_PyDiffrecObj * const self
)
{
    if (!self)
        return;

    if (!self->tau)
        return;

    if (self->tau->error_occurred)
        return;

    if (self->verbose)
        puts("\r\tDiffractionCorrection: Passing keywords to C struct...");

    self->tau->sigma = self->sigma;
    self->tau->bfac = self->bfac;
    self->tau->eccentricity = self->eccentricity;
    self->tau->periapse = self->periapse;
    self->tau->use_fwd = self->use_fwd;
    self->tau->use_norm = self->use_norm;
}
