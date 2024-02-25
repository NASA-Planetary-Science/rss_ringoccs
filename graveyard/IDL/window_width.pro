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
;       Computes the window width for a given resolution and geometry.         ;
;   Arguments:                                                                 ;
;       RES (real):                                                            ;
;           The requested resolution.                                          ;
;       NORMEQ (real):                                                         ;
;           The normalized equivalent width of the desired window function.    ;
;       FSKY (real, array-like):                                               ;
;           The sky frequency.                                                 ;
;       FRES (real, array-like):                                               ;
;           The Fresnel scale.                                                 ;
;       RHO_DOT (real, array-like):                                            ;
;           The time derivative of the radius of the point in the plane.       ;
;   Keywords:                                                                  ;
;       SIGMA (real):                                                          ;
;           The Allen deviation.                                               ;
;   Output:                                                                    ;
;       W_KM_VALS (real, array-like):                                          ;
;           The window width for each of the points.                           ;
;   References:                                                                ;
;       1.) Marouf, Tyler, Rosen, 1986                                         ;
;           Profiling Saturn's Rings by Radio Occultation.                     ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Computes the window width for a given resolution, window type, and geometry.
FUNCTION WINDOW_WIDTH, RES, NORMEQ, FSKY, FRES, RHO_DOT, SIGMA = SIGMA

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; If no user-requested Allen deviation is provided use Cassini's.
    IF KEYWORD_SET(SIGMA) THEN BEGIN
        SIGMA = SIGMA
    ENDIF ELSE BEGIN
        SIGMA = 2.0E-13
    ENDELSE

    ; Parameters for the Lambert W function and window window.
    OMEGA = 2.0 * !PI * FSKY
    ALPHA_FACTOR = OMEGA * SIGMA
    ALPHA = ALPHA_FACTOR * ALPHA_FACTOR / (2.0 * RHO_DOT)
    RCPR_ALPHA = 1.0 / ALPHA
    P = RES * RCPR_ALPHA / (FRES * FRES)

    ; The inverse exists only if P > 1.
    IF MIN(P) LT 1.0001 THEN BEGIN
        MESSAGE, 'Illegal parameters. Exclude points or use coarser resolution.'
    ENDIF

    ; The inverse formula for resolution-to-window-width is given in terms of
    ; the Lambert W function. This is computed via INVERSE_RESOLUTION.
    W_KM_VALS = SQRT(2.0) * NORMEQ * INVERSE_RESOLUTION(P) * RCPR_ALPHA
    RETURN, W_KM_VALS
END
