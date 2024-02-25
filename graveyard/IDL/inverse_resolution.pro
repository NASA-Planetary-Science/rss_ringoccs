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
;       Computes the inverse of y = x / (exp(-x) + x - 1) for y > 1.           ;
;   Arguments:                                                                 ;
;       XIN (real, array-like):                                                ;
;           The independent variable.                                          ;
;   Output:                                                                    ;
;       INV_RES (real, array-like):                                            ;
;           The inverse of y = x / (exp(-x) + x - 1).                          ;
;   References:                                                                ;
;       1.) Marouf, Tyler, Rosen, 1986                                         ;
;           Profiling Saturn's Rings by Radio Occultation.                     ;
;   Notes:                                                                     ;
;       The function occurs in the equation for the resolution of the          ;
;       diffraction reconstruction when the Allen deviation is taken into      ;
;       account, hence the name.                                               ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Function for computing the inverse of x / (exp(-x) + x - 1).
FUNCTION INVERSE_RESOLUTION, XIN

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; The inverse uses Exp and Lambert W. Compute these terms.
    P1 = XIN / (1.0 - XIN)
    P2 = P1 * EXP(P1)

    ; The final expression can be written in terms of the Lambert W function.
    RETURN, LAMBERTW(P2) - P1
END
