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
;       Numerically computes the second derivative of a function.              ;
;   Arguments:                                                                 ;
;       FUNC (real, array-like):                                               ;
;           A 5xN array representing the stencil data for a function.          ;
;       DX (REAL):                                                             ;
;           The displacement between samples in the five-point stencil.        ;
;   Output:                                                                    ;
;       DOUBLE_PRIME (real, array-like):                                       ;
;           The second derivative of FUNC.                                     ;
;   References:                                                                ;
;       1.) Abramowitz, Stegun, 1970                                           ;
;           Handbook of Mathematical Functions.                                ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Function for approximating second derivates using a five-point stencil.
FUNCTION LAGRANGE_D2, FUNC, DX

    ; Error checking codes.
    ON_ERROR, 2

    ; Info about the data in the func variable.
    SIZE_FUNC = SIZE(F)

    ; The input should be 5xN where N is any positive integer.
    IF (SIZE_FUNC[0] NE 2) OR (SIZE_FUNC[1] NE 5) THEN BEGIN
        MESSAGE, 'Invalid Input. FUNC Must be a 5xN Array'
    ENDIF

    ; Denominator of the Lagrange formula.
    DEN = 12.0 * DX * DX

    ; The numerator, weighted sum of the stencil points.
    NUM = -FUNC[0, *]+16.0*FUNC[1, *]-30.0*FUNC[2, *]+16.0*FUNC[3, *]-FUNC[4, *]

    ; The approximation for the second derivative.
    F_DOUBLE_PRIME = NUM / DEN
    RETURN, F_DOUBLE_PRIME
END
