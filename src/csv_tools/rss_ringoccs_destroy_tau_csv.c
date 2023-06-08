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
 *      Function for free'ing a Tau CSV object and free'ing all of the        *
 *      pointers contained inside the struct.                                 *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 31, 2020                                             *
 ******************************************************************************/

/*  free is found here, as is NULL.                                           */
#include <stdlib.h>

/*  rssringoccs_TauCSV typedef here, and function prototype given.            */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Function for freeing the memory of a TauCSV object.                       */
void rssringoccs_Destroy_TauCSV(rssringoccs_TauCSV **tau)
{
    /*  Variable for a pointer to the TauCSV object.                          */
    rssringoccs_TauCSV *tau_inst;

    /*  If the input pointer is NULL, there is nothing to free.               */
    if (tau == NULL)
        return;

    /*  Otherwise, get a pointer to the TauCSV object.                        */
    tau_inst = *tau;

    /*  If this pointer is NULL, there is nothing to free.                    */
    if (tau_inst == NULL)
        return;

    /*  Free all of the members of the TauCSV object.                         */
    rssringoccs_Destroy_TauCSV_Members(tau_inst);

    /*  If an error occurred, the error_message variable is malloced and a    *
     *  string is stored in it. Free this pointer if this is the case.        */
    if (tau_inst->error_message != NULL)
    {
        free(tau_inst->error_message);

        /*  Set the pointer to NULL to prevent freeing it twice.              */
        tau_inst->error_message = NULL;
    }

    /*  Free the TauCSV pointer and set to NULL to prevent freeing it twice.  */
    free(tau_inst);
    *tau = NULL;
    return;
}
/*  End of rssringoccs_Destroy_TauCSV.                                        */
