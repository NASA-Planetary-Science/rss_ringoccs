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
 *      Function for free'ing a Cal CSV object and free'ing all of the        *
 *      pointers contained inside the struct.                                 *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 31, 2020                                             *
 ******************************************************************************/

/*  Macro for freeing a pointer and setting it to NULL.                       */
#include <libtmpl/include/compat/tmpl_free.h>

/*  rssringoccs_CalCSV typedef provided here.                                 */
#include <rss_ringoccs/include/types/rss_ringoccs_calcsv.h>

/*  NULL macro defined here.                                                  */
#include <stddef.h>

/*  Function prototype / forward declaration.                                 */
extern void rssringoccs_CalCSV_Destroy(rssringoccs_CalCSV ** const cal);

/*  Tell the compiler about the main destructor function.                     */
extern void rssringoccs_CalCSV_Destroy_Members(rssringoccs_CalCSV * const cal);

/*  Function for freeing the memory in a CalCSV object.                       */
void rssringoccs_CalCSV_Destroy(rssringoccs_CalCSV ** const cal)
{
    /*  Used for the pointer to the CSV object.                               */
    rssringoccs_CalCSV *cal_inst;

    /*  If the input pointer is NULL, simply return.                          */
    if (!cal)
        return;

    /*  Otherwise, get a pointer to the CalCSV object.                        */
    cal_inst = *cal;

    /*  If this is NULL, there's no need to free it. Return.                  */
    if (!cal_inst)
        return;

    /*  Free all of the pointers inside the CalCSV object.                    */
    rssringoccs_CalCSV_Destroy_Members(cal_inst);

    /*  The error_message member is a pointer to a constant string, it does   *
     *  not need to be freed. Set the pointer to NULL to avoid reading it.    */
    cal_inst->error_message = NULL;

    /*  Free the CalCSV pointer and set it to NULL to prevent freeing twice.  */
    TMPL_FREE(*cal);
}
/*  End of rssringoccs_CalCSV_Destroy.                                        */
