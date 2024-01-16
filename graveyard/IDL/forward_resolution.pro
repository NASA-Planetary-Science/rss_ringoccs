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
;       Computes the resolution factor from the Allen deviation.               ;
;   Arguments:                                                                 ;
;       ALDEV (real, array-like):                                              ;
;           The Allen deviation.                                               ;
;   Output:                                                                    ;
;       RES (real, array-like):                                                ;
;           The resolution factor for the Allen deviation.                     ;
;   References:                                                                ;
;       1.) Marouf, Tyler, Rosen, 1986                                         ;
;           Profiling Saturn's Rings by Radio Occultation.                     ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Function for computing the inverse of x / (exp(-x) + x - 1).
FUNCTION FORWARD_RESOLUTION, ALDEV

    ; Error checking code.
    ON_ERROR, 2

    ; Ratio of two terms. Compute the numerator and denominator.
    NUM = 0.5 * ALDEV * ALDEV
    DEN = EXP(-ALDEV) + ALDEV - 1.0

    ; The resolution factor is the ratio of these terms.
    RETURN, NUM / DEN
END
