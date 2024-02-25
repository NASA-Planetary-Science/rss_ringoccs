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
;       Computes the Kaiser-Bessel window with alpha = 3.5 pi.                 ;
;   Arguments:                                                                 ;
;       WIDTH (real):                                                          ;
;           The window width.                                                  ;
;       DX (real):                                                             ;
;           The spacing between samples.                                       ;
;   Output:                                                                    ;
;       W_FUNC (real, array-like):                                             ;
;           The window function with evenly spaced samples.                    ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Computes the Kaiser-Bessel window function with alpha = 3.5 pi.
FUNCTION WINDOW_ARRAY_KB35, WIDTH, DX

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; The width should be a single positive number.
    IF (CHECK_POSITIVE_REAL_NUMBER(WIDTH) EQ 0) THEN BEGIN
        MESSAGE, "Window width must be a single positive real number."
    ENDIF

    ; The displacement between samples should also be a single positive number.
    IF (CHECK_POSITIVE_REAL_NUMBER(DX) EQ 0) THEN BEGIN
        MESSAGE, "Sample width must be a single positive real number."
    ENDIF

    ; The number of points in the window. Should be an odd number.
    NW_PTS = 2*FLOOR(WIDTH / (2.0*DX)) + 1

    ; The independent variable, evenly spaced samples.
    XARR = (FINDGEN(NW_PTS) - (NW_PTS - 1)*0.5) * DX

    ; The scale factors for the argument and window.
    ALPHA_PI = 3.5 * !PI
    BESSEL_ALPHA_PI = 7257.8007901476148975659063117440177659

    ; The argument for the window function, and the factor for the sqrt term.
    ARG = 2.0 * XARR / WIDTH
    FACTOR = 1.0 - ARG*ARG

    ; The formula for the Kaiser-Bessel window.
    W_FUNC = BESELI(ALPHA_PI * SQRT(FACTOR), 0.0) / BESSEL_ALPHA_PI
    RETURN, W_FUNC
END
