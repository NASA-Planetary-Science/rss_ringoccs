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
;       Computes the cosine-squared window function.                           ;
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

; Computes the cosine-squared window function.
FUNCTION WINDOW_ARRAY_COSS, WIDTH, DX

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

    ; The argument for the cosine function.
    ARG = !PI * XARR / WIDTH

    ; The formula for the cosine-squared window.
    COS_ARG = COS(ARG)
    W_FUNC = COS_ARG * COS_ARG
    RETURN, W_FUNC
END
