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
 *      Function for free'ing a DLP CSV object and free'ing all of the        *
 *      pointers contained inside the struct.                                 *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 31, 2020                                             *
 ******************************************************************************/

/*  free is found here, as is NULL.                                           */
#include <stdlib.h>

/*  rssringoccs_DLPCSV typedef here, and function prototype given.            */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Function for freeing the memory in a DLPCSV object.                       */
void rssringoccs_Destroy_DLPCSV(rssringoccs_DLPCSV **dlp)
{
    /*  Variable for a pointer to the DLPCSV object.                          */
    rssringoccs_DLPCSV *dlp_inst;

    /*  If the input pointer is NULL, do not attempt to free it.              */
    if (dlp == NULL)
        return;

    /*  Get a pointer to the DLPCSV object.                                   */
    dlp_inst = *dlp;

    /*  If this pointer is NULL, there is nothing to free.                    */
    if (dlp_inst == NULL)
        return;

    /*  Free all of the members inside the DLP object.                        */
    rssringoccs_Destroy_DLPCSV_Members(dlp_inst);

    /*  If an error occurred, the error_message variable is malloced and a    *
     *  string is stored inside of it. Free this pointer if this is true.     */
    if (dlp_inst->error_message != NULL)
    {
        free(dlp_inst->error_message);

        /*  Set the pointer to NULL to avoid freeing it twice.                */
        dlp_inst->error_message = NULL;
    }

    /*  Free the DLPCSV object pointer and set it to NULL to prevent trying   *
     *  to free the pointer twice.                                            */
    free(dlp_inst);
    *dlp = NULL;
    return;
}
/*  End of rssringoccs_Destroy_DLPCSV.                                        */
