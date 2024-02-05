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
;       Computes the modified Kaiser-Bessel window with variable alpha.        ;
;   Arguments:                                                                 ;
;       WIDTH (real):                                                          ;
;           The window width.                                                  ;
;       DX (real):                                                             ;
;           The spacing between samples.                                       ;
;       ALPHA (real):                                                          ;
;           The alpha parameter for the Kaiser-Bessel function.                ;
;   Output:                                                                    ;
;       W_FUNC (real, array-like):                                             ;
;           The window function with evenly spaced samples.                    ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Computes the modified Kaiser-Bessel window function with positive alpha.
FUNCTION WINDOW_ARRAY_KBMD, WIDTH, DX, ALPHA

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

    ; Same check for the alpha parameter.
    IF (CHECK_POSITIVE_REAL_NUMBER(ALPHA) EQ 0) THEN BEGIN
        MESSAGE, "Alpha factor must be a single positive real number."
    ENDIF

    ; The number of points in the window. Should be an odd number.
    NW_PTS = 2*FLOOR(WIDTH / (2.0*DX)) + 1

    ; The independent variable, evenly spaced samples.
    XARR = (FINDGEN(NW_PTS) - (NW_PTS - 1)*0.5) * DX

    ; The scale factors for the argument and window.
    ALPHA_PI = ALPHA * !PI
    BESSEL_ALPHA_PI = BESELI(ALPHA_PI, 0.0)

    ; The argument for the window function, and the factor for the sqrt term.
    ARG = 2.0 * XARR / WIDTH
    FACTOR = 1.0 - ARG*ARG

    ; Numerator and denominator are modified from the standard KB function.
    NUM = BESELI(ALPHA_PI * SQRT(FACTOR), 0.0) - 1.0
    DEN = BESSEL_ALPHA_PI - 1.0

    ; The formula for the modified Kaiser-Bessel window is the ratio of these.
    W_FUNC = NUM / DEN
    RETURN, W_FUNC
END
