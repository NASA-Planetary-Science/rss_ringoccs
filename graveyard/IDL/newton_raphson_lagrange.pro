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
;       Computes the Newton iterate of the derivative of a function.           ;
;   Arguments:                                                                 ;
;       FUNC (real, array-like):                                               ;
;           A 5xN array representing the stencil data for a function.          ;
;       DX (real):                                                             ;
;           The displacement between samples in the five-point stencil.        ;
;   Output:                                                                    ;
;       NEWTON_ITER (real, array-like):                                        ;
;           The value -f' / f'' of the input function f.                       ;
;   References:                                                                ;
;       1.) Abramowitz, Stegun, 1970                                           ;
;           Handbook of Mathematical Functions.                                ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018/04/18                                                         ;
;------------------------------------------------------------------------------;

; Numerically computes the Newton-Iterate of the derivative of a function.
FUNCTION NEWTON_RAPHSON_LAGRANGE, FUNC, DX

    ; Error checking code.
    ON_ERROR, 2

    ; Info about the data in the func variable.
    SIZE_FUNC = SIZE(FUNC)

    ; The input should be 5xN where N is any positive integer.
    IF (SIZE_FUNC[0] NE 2) OR (SIZE_FUNC[1] NE 5) THEN BEGIN
        MESSAGE, 'Invalid Input. FUNC must be a 5xN Array'
    ENDIF

    ; Use Lagrange's method to compute the first derivative.
    F_PRIME = LAGRANGE_D1(FUNC, DX)

    ; Similar idea for the numerical second derivative.
    F_DOUBLE_PRIME = LAGRANGE_D2(FUNC, DX)

    ; Newton iterate, -f' / f''. Compute this.
    ROOT = -F_PRIME / F_DOUBLE_PRIME
    RETURN, ROOT
END
