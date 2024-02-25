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
;       Creates a triangular window for a data set.                            ;
;   Arguments:                                                                 ;
;       XIN (real, array-like):                                                ;
;           The independent variable.                                          ;
;   Output:                                                                    ;
;       TRIANGLE (real, array-like):                                           ;
;           A triangular window for XIN.                                       ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Function for creating a triangular window for a data set.
FUNCTION TRIANGULAR_WINDOW, XIN

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; The input should be real-valued.
    IF (CHECK_REAL(XIN) EQ 0) THEN BEGIN
        MESSAGE, "Input should be a real-valued array."
    ENDIF

    ; Size data for the input argument.
    X_LENGTH = N_ELEMENTS(XIN)
    HALF_X_LENGTH = X_LENGTH / 2

    ; Factor for the triangle function.
    FACTOR = 2.0 / FLOAT(X_LENGTH)

    ; Create a straight line at first.
    TRIANGLE = FINDGEN(X_LENGTH) * FACTOR

    ; Then flip the second half of the line to make a triangle.
    TRIANGLE[HALF_X_LENGTH:-1] = 2.0 - TRIANGLE[HALF_X_LENGTH:-1]
    RETURN, TRIANGLE
END
