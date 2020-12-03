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
 *                        rss_ringoccs_idl_bessel_j0                          *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Provide an example of using the C-IDL API to call librssringoccs      *
 *      functions in IDL.                                                     *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_special_functions.h:                                     *
 *          Library of special functions written in C.                        *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 3, 2020                                              *
 ******************************************************************************/

#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

/*  The function will be callable from IDL via the string                     *
 *  'rssringoccs_IDL_Bessel_J0' using the CALL_EXTERNAL routine.              */
void rssringoccs_IDL_Bessel_J0(int argc, void *argv[])
{
    /*  Declare two double pointers, the input and output for IDL.            */
    double *x, *y;

    /*  Size is the size of the input IDL array.                              */
    long n, size;

    /*  Get the parameters passed from IDL.                                   */
    x    = (double *)argv[0];
    size = *(long *) argv[1];
    y    = (double *)argv[2];

    /*  Loop through each point and compute the Bessel J0 function.           */
    for (n=0; n<size; ++n)
        y[n] = rssringoccs_Double_Bessel_J0(x[n]);
}

/******************************************************************************
 *  Compile this with:                                                        *
 *      gcc rss_ringoccs_idl_bessel_j0.c -O3 -shared -o                       *
 *          rss_ringoccs_idl_bessel_j0.so -lrssringoccs                       *
 *                                                                            *
 *  -shared means we're creating a shared object (library file).              *
 *  -O3 means we want level 3 optimization to make it faster.                 *
 *  -o rss_ringoccs_idl_bessel_j0.so means we're compiling this into a file   *
 *      called 'rss_ringoccs_idl_bessel_j0.so'.                               *
 *  -lrssringoccs is a linker parameter, telling the linking to link the      *
 *      program to librssringoccs.so which should be in /usr/local/lib/.      *
 *                                                                            *
 *  See bessel_j0.pro for details on the IDL part of this program.            *
 ******************************************************************************/
