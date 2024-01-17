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
;       Numerically computes the Newton iterator of a functions derivative.    ;
;   Arguments:                                                                 ;
;       FUNC (string):                                                         ;
;           The function name.                                                 ;
;       XIN (real, complex, array-like):                                       ;
;           The independent variable for the function.                         ;
;       DX (REAL):                                                             ;
;           The displacement between samples in the five-point stencil.        ;
;   Output:                                                                    ;
;       NEWTON_ITER (real, complex, array-like):                               ;
;           The value -f' / f''.                                               ;
;   References:                                                                ;
;       1.) Abramowitz, Stegun, 1970                                           ;
;           Handbook of Mathematical Functions.                                ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018                                                               ;
;------------------------------------------------------------------------------;

; Numerical Newton's method using five-point stencils.
FUNCTION NEWTON_RAPHSON_LAGRANGE_FROM_FUNC, FUNC, XIN, DX

    ; Error checking code.
    ON_ERROR, 2

    ; The input function should be a string.
    IF (IDL_TYPE(FUNC) NE IDL_TYPE("String")) THEN BEGIN
        MESSAGE, "Input function should be a string (the name of the function)."
    ENDIF

    ; The data should be real or complex.
    IF (CHECK_REAL(XIN) EQ 0) AND (CHECK_COMPLEX(XIN) EQ 0) THEN BEGIN
        MESSAGE, "Input array should be real or complex valued."
    END

    ; The displacement should be strictly positive.
    IF (CHECK_POSITIVE_REAL_NUMBER(DX) EQ 0) THEN BEGIN
        MESSAGE, "Displacement should be a strictly positive real number"
    END

    ; The stencil points. The third is just XIN.
    X0 = XIN - 2.0 * DX
    X1 = XIN - 1.0 * DX
    X2 = XIN
    X3 = XIN + 1.0 * DX
    X4 = XIN + 2.0 * DX

    ; The function evaluated at the stencil points.
    Y0 = CALL_FUNCTION(FUNC, X0)
    Y1 = CALL_FUNCTION(FUNC, X1)
    Y2 = CALL_FUNCTION(FUNC, X2)
    Y3 = CALL_FUNCTION(FUNC, X3)
    Y4 = CALL_FUNCTION(FUNC, X4)

    ; Denominators of the Lagrange formulae.
    DEN1 = 12.0 * DX
    DEN2 = DEN1 * DX

    ; The numerators, weighted sums of the stencil points.
    NUM1 = Y0 - 8.0*Y1 + 8.0*Y3 - Y4
    NUM2 = -Y0 + 16.0*Y1 - 30.0*Y2 + 16.0*Y3 - Y4

    ; The first and second derivative.
    F_PRIME = NUM1 / DEN1
    F_DOUBLE_PRIME = NUM2 / DEN2

    ; The Newton iterate for the derivative of f.
    NEWTON_ITER = -F_PRIME / F_DOUBLE_PRIME
    RETURN, NEWTON_ITER
END
