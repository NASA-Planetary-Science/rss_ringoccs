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
 ******************************************************************************/
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Standard algorithm of time-complexity O(n) to take an array and reverse   *
 *  the order. This is equivalent to taking a numpy array arr and writing     *
 *  arr = arr[::-1] (in Python 3, at least). Since all of the main routines   *
 *  use double pointers, only a double version is provided.                   */
void rssringoccs_Reverse_Float_Array(float *arr, unsigned long arrsize)
{
    float val;
    unsigned long i;
    for(i=0; i<arrsize/2; i++)
    {
        val = arr[i];
        arr[i] = arr[arrsize-i-1];
        arr[arrsize-i-1] = val;
    }
}

void rssringoccs_Reverse_Double_Array(double *arr, unsigned long arrsize)
{
    double val;
    unsigned long i;
    for(i=0; i<arrsize/2; i++)
    {
        val = arr[i];
        arr[i] = arr[arrsize-i-1];
        arr[arrsize-i-1] = val;
    }
}

void rssringoccs_Reverse_LDouble_Array(long double *arr, unsigned long arrsize)
{
    long double val;
    unsigned long i;
    for(i=0; i<arrsize/2; i++)
    {
        val = arr[i];
        arr[i] = arr[arrsize-i-1];
        arr[arrsize-i-1] = val;
    }
}
