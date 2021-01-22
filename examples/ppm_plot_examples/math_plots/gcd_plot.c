/******************************************************************************
 *                                 LICENSE                                    *
 ******************************************************************************
 *  This file is free software: you can redistribute it and/or modify it      *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  This is distributed in the hope that it will be useful,                   *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 ******************************************************************************/

#include <stdio.h>

/*  Simple routine for computing the GCD of non-negative numbers.             */
long GCD(long n0, long n1)
{
    /*  Special case. The while loop can be infinite if one of the entries is *
     *  zero. GCD(n, 0) = n, so use this.                                     */
    if (n0 == 0)
        return n1;
    else if (n1 == 0)
        return n0;

    while(n0 != n1)
    {
        if(n0 > n1)
            n0 -= n1;
        else
            n1 -= n0;
    }
    return n0;
}
/*  End of GCD.                                                               */

/*  Function for computing the "area" of the blocks with GCD(x, y) = 1. This  *
 *  simply sums over the square [0, n] x [0, n] for which entries have GCD 1, *
 *  and then divides by N^2.                                                  */
double test(unsigned int N)
{
    unsigned int x, y;
    unsigned long counter = 0UL;

    for (x=0; x<N; ++x)
        for (y=0; y<N; ++y)
            if (GCD(x, y) == 1)
                ++counter;

    return (double)counter / (double)(N*N);
}
/*  End of "test" function.                                                   */

/*  Function for plotting GCD grid. White for GCD(x, y) = 1, black otherwise. */
int main(void)
{
    /*  Declare all necessary variables.                                      */
    unsigned int x, y, z_x, z_y;
    unsigned N = 1024U;
    unsigned int n;
    unsigned char black = 0U;
    unsigned char white = 255U;

    /*  I only want to draw the block [0, 63] x [0, 63], but a 64x64 PGM file *
     *  will be really-really small and zooming in makes it blurry. Use this  *
     *  scale factor to draw the [0, 63] x [0, 63] region in a 1024x1024 PGM. */
    double scale = 64.0 / N;

    /*  And a variable for writing to a file.                                 */
    FILE *fp, *ftxt;

    /*  Open the PGM file and write the preamble to it.                       */
    fp = fopen("gcd_plot.pgm", "w");
    fprintf(fp, "P5\n%d %d\n255\n", N, N);

    /*  Loop through each pixel.                                              */
    for (y=0; y<N; ++y)
    {
        /*  PGMs plot top-to-bottom, whereas mathematicians think in a        *
         *  bottom-to-top style. Plot with the z_y variable to fix this.      */
        z_y = (unsigned int)(scale*(N - y));

        /*  Loop over every x pixel as well.                                  */
        for (x=0; x<N; ++x)
        {
            z_x = (unsigned int)(scale*x);

            /*  If z_x and z_y are coprime, plot the pixel white.             */
            if (GCD(z_x, z_y) == 1)
                fputc(white, fp);
            else
                fputc(black, fp);
        }
        /*  End of for-loop for the pixels x-coordinate.                      */
    }
    /*  End of for-loop for the pixels y-coordinate.                          */

    /*  Close the file.                                                       */
    fclose(fp);

    /*  Run the test function on various sizes.                               */
    N = 1000;

    /*  We'll write the output to a text file so we can make neater plots     *
     *  using either GNU plotutils or matplotlib in Python.                   */
    ftxt = fopen("gcd_test.txt", "w");

    /*  n = 0 is rather pointless, so start at n = 1.                         */
    for (n=1; n<=N; ++n)
        fprintf(fp, "%.16f\n", test(n));

    /*  Close the file.                                                       */
    fclose(ftxt);

    /*  And lastly, let's test a bigger value.                                */
    printf("%.16f\n", test(2E4));
    return 0;
}
/*  End of main.                                                              */
