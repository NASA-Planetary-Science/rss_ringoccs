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
;       Computes the numerator of the normalized equivalent width.             ;
;   Arguments:                                                                 ;
;       WINDOW_ARRAY (real, array-like):                                       ;
;           A data set corresponding to the values of the window function.     ;
;   Output:                                                                    ;
;       WSQ_HALF (real):                                                       ;
;           The numerator of the normalized equivalent width.                  ;
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

; Computes the numerator of the normalized equivalent width.
FUNCTION CALC_WSQ_HALF, WINDOW_ARRAY

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Only positive values in the array are needed.
    WINDOW_ARRAY_POSITIVE_INDICES = WHERE(WINDOW_ARRAY GT 0.d0)

    ; Subarray of the window array with the positive values.
    WINDOW_ARRAY_POSITIVE = WINDOW_ARRAY[WINDOW_ARRAY_POSITIVE_INDICES]

    ; The numerator of the normalized equivalent width.
    RETURN, MEAN(WINDOW_ARRAY_POSITIVE^2)
END
