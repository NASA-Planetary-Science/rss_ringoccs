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
 *      Function for free'ing a Geo CSV object and free'ing all of the        *
 *      pointers contained inside the struct.                                 *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 31, 2020                                             *
 ******************************************************************************/

/*  free is found here, as is NULL.                                           */
#include <stdlib.h>

/*  rssringoccs_GeoCSV typedef here, and function prototype given.            */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Function for freeing all of the memory in a GeoCSV struct.                */
void rssringoccs_GeoCSV_Destroy(rssringoccs_GeoCSV **geo)
{
    /*  Variable for a pointer to the GeoCSV object.                          */
    rssringoccs_GeoCSV *geo_inst;

    /*  If the input pointer is NULL, there's no need to free it.             */
    if (!geo)
        return;

    /*  Otherwise, get a pointer to the GeoCSV object.                        */
    geo_inst = *geo;

    /*  If this pointer is NULL, there's no need to free anything.            */
    if (!geo_inst)
        return;

    /*  Free all of the members in the GeoCSV object.                         */
    rssringoccs_GeoCSV_Destroy_Members(geo_inst);

    /*  If an error occurred, error_message is malloced and a string is       *
     *  stored in the variable. Free the pointer if this is the case.         */
    if (geo_inst->error_message != NULL)
    {
        free(geo_inst->error_message);

        /*  Set the pointer to NULL to prevent freeing it twice.              */
        geo_inst->error_message = NULL;
    }

    /*  Free the GeoCSV pointer and set it to NULL to prevent freeing twice.  */
    free(geo_inst);
    *geo = NULL;
    return;
}
/*  End of rssringoccs_GeoCSV_Destroy.                                        */
