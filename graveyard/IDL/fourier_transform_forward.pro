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
;       Computes the lazy O(n^2) forward Fourier transform.                    ;
;   Arguments:                                                                 ;
;       X_IN (real, array-like):                                               ;
;           The independent variable.                                          ;
;       F (complex, array-like):                                               ;
;           The data set.                                                      ;
;   Output:                                                                    ;
;       F_HAT (complex, array-like):                                           ;
;           The transformed data.                                              ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Function for computing the Fourier transform of a data set.
FUNCTION FOURIER_TRANSFORM_FORWARD, X_IN, F

    ; Step size between samples. Assumed uniform.
    DX = X_IN[1] - X_IN[0]

    ; Matrix product from the input, represents xt with transform variable t.
    X_IN_SQ = X_IN # X_IN

    ; The Fourier kernel in the transform.
    OMEGA = COMPLEX(0.0, -2.0 * !PI * X_IN_SQ)
    KERNEL = EXP(OMEGA)

    ; The discrete transform can be computed from a matrix product. This is
    ; the lazy method and runs in O(n^2) instead of O(n log(n)).
    F_HAT = KERNEL # F * DX
    RETURN, F_HAT
END
