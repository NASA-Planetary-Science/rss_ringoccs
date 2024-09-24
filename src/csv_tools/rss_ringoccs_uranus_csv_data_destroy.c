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
 *      Function for free'ing a CSV object and free'ing all of the            *
 *      pointers contained inside the struct.                                 *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 31, 2020                                             *
 ******************************************************************************/

/*  free is found here, as is NULL.                                           */
#include <stdlib.h>

/*  rssringoccs_CSVData typedef here, and function prototype given.           */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Function for freeing the memory in a CSV object.                          */
void rssringoccs_UranusCSVData_Destroy(rssringoccs_UranusCSVData **csv)
{
    /*  Used for the pointer to the CSV object.                               */
    rssringoccs_UranusCSVData *csv_inst;

    /*  If the input pointer is NULL, simply return.                          */
    if (csv == NULL)
        return;

    /*  Otherwise, get a pointer to the CSV object.                           */
    csv_inst = *csv;

    /*  If this is NULL, there's no need to free it. Return.                  */
    if (csv_inst == NULL)
        return;

    /*  Free all of the pointers inside the CSV object.                       */
    rssringoccs_UranusCSVData_Destroy_Members(csv_inst);

    /*  If an error occurred along the way, the error_message variable is     *
     *  malloced and a string is stored. Check if we need to free this.       */
    if (csv_inst->error_message != NULL)
    {
        free(csv_inst->error_message);

        /*  To avoid freeing twice, reset the pointer to NULL.                */
        csv_inst->error_message = NULL;
    }

    /*  Free the pointer to the object and set it to NULL to avoid freeing    *
     *  this object twice.                                                    */
    free(csv_inst);
    *csv = NULL;
    return;
}
/*  End of rssringoccs_UranusCSVData_Destroy.                                 */
