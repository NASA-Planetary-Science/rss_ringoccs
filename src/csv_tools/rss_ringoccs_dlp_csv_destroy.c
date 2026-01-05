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

/*  Macro for freeing a pointer and setting it to NULL.                       */
#include <libtmpl/include/compat/tmpl_free.h>

/*  rssringoccs_DLPCSV typedef provided here.                                 */
#include <rss_ringoccs/include/types/rss_ringoccs_dlpcsv.h>

/*  NULL macro defined here.                                                  */
#include <stddef.h>

/*  Function prototype / forward declaration.                                 */
extern void rssringoccs_DLPCSV_Destroy(rssringoccs_DLPCSV ** const dlp);

/*  Tell the compiler about the main destructor function.                     */
extern void rssringoccs_DLPCSV_Destroy_Members(rssringoccs_DLPCSV * const dlp);

/*  Function for freeing the memory in a DLPCSV object.                       */
void rssringoccs_DLPCSV_Destroy(rssringoccs_DLPCSV **dlp)
{
    /*  Variable for a pointer to the DLPCSV object.                          */
    rssringoccs_DLPCSV *dlp_inst;

    /*  If the input pointer is NULL, do not attempt to free it.              */
    if (!dlp)
        return;

    /*  Get a pointer to the DLPCSV object.                                   */
    dlp_inst = *dlp;

    /*  If this pointer is NULL, there is nothing to free.                    */
    if (!dlp_inst)
        return;

    /*  Free all of the members inside the DLP object.                        */
    rssringoccs_DLPCSV_Destroy_Members(dlp_inst);

    /*  The error_message member is a pointer to a constant string, it does   *
     *  not need to be freed. Set the pointer to NULL to avoid reading it.    */
    dlp_inst->error_message = NULL;

    /*  Free the DLPCSV object pointer and set it to NULL to prevent trying   *
     *  to free the pointer twice.                                            */
    TMPL_FREE(*dlp);
}
/*  End of rssringoccs_DLPCSV_Destroy.                                        */
