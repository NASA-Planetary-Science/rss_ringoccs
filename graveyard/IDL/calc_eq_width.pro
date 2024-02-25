;------------------------------------------------------------------------------;
;                                   LICENSE                                    ;
;------------------------------------------------------------------------------;
;   This file is part of rss_ringoccs.                                         ;
;                                                                              ;
;   rss_ringoccs is free software: you can redistribute it and/or              ;
;   modify it under the terms of the GNU General Public License as published   ;
;   by the Free Software Foundation, either version 3 of the License, or       ;
;   (at your option) any later version.                                        ;
;                                                                              ;
;   rss_ringoccs is distributed in the hope that it will be useful             ;
;   but WITHOUT ANY WARRANTY; without even the implied warranty of             ;
;   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              ;
;   GNU General Public License for more details.                               ;
;                                                                              ;
;   You should have received a copy of the GNU General Public License          ;
;   along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.     ;
;------------------------------------------------------------------------------;
;   Purpose:                                                                   ;
;       Computes the normalized equivalent width of a window function.         ;
;   Arguments:                                                                 ;
;       WINDOW_ARRAY (real, array-like):                                       ;
;           A data set corresponding to the values of the window function.     ;
;   Output:                                                                    ;
;       EQ_WIDTH (real):                                                       ;
;           The normalized equivalent width.                                   ;
;   References:                                                                ;
;       1.) Marouf, Tyler, Rosen, 1986                                         ;
;           Profiling Saturn's Rings by Radio Occultation.                     ;
;       2.) Carroll, Ostie, 2007                                               ;
;           An Introduction to Modern Astrophysics.                            ;
;   Notes:                                                                     ;
;       1.) See equation 20 of reference (1) for more details.                 ;
;       2.) It is assumed your window array has x-values between -0.5 and 0.5. ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Calculate the normalized equivalent width of a windowing function.
FUNCTION CALC_EQ_WIDTH, WINDOW_ARRAY

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Only the positive values of the array contribute to the integrals.
    WINDOW_ARRAY_POSITIVE_INDICES = WHERE(WINDOW_ARRAY GT 0.d0)

    ; Subarray with only positive values.
    WINDOW_ARRAY_POSITIVE = WINDOW_ARRAY[WINDOW_ARRAY_POSITIVE_INDICES]

    ; The mean of the square of the array, mean(w^2).
    WINDOW_MEAN_SQUARE = MEAN(WINDOW_ARRAY_POSITIVE^2)

    ; The square of the mean of the array, mean(w)^2.
    WINDOW_SQUARE_MEAN = MEAN(WINDOW_ARRAY_POSITIVE)^2

    ; The equivalent width of the window function is the ratio of these two.
    RETURN, WINDOW_MEAN_SQUARE / WINDOW_SQUARE_MEAN
END
