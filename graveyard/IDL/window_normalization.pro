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
;       Computes the normalizing factor of a given window function.            ;
;   Arguments:                                                                 ;
;       RHO (real, array-like):                                                ;
;           The radial data.                                                   ;
;       W_FUNC (real, array-like):                                             ;
;           The window data.                                                   ;
;       F_SCALE (real):                                                        ;
;           The Fresnel scale.                                                 ;
;   Output:                                                                    ;
;       NORM_FACTOR (real):                                                    ;
;           The normalization factor.                                          ;
;   Notes:                                                                     ;
;       The Fresnel transform is given by the integral of T * w * exp(i psi)   ;
;       where w is the window function, psi is the Fresnel kernel, and T is    ;
;       the actual data. The window function makes free-space regions have     ;
;       transforms that are less than 1 in magnitude. By dividing by           ;
;       the integral of w * exp(i psi), which is the free space integral, we   ;
;       ensure free space has a transmittance with magnitude 1.                ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018                                                               ;
;------------------------------------------------------------------------------;
FUNCTION WINDOW_NORMALIZATION, RHO, W_FUNC, F_SCALE

    ; Error checking code.
    ON_ERROR, 2

    ; The radial data must be real.
    IF (CHECK_REAL(RHO) EQ 0) THEN BEGIN
        MESSAGE, "RHO must be real valued"
    ENDIF

    ; So must the window function.
    IF (CHECK_REAL(W_FUNC) EQ 0) THEN BEGIN
        MESSAGE, "W_FUNC must be real valued"
    ENDIF

    ; The Fresnel scale is also real-valued.
    IF (CHECK_REAL(F_SCALE) EQ 0) THEN BEGIN
        MESSAGE, "F_SCALE must be a positive real number"
    ENDIF

    ; The number of points in the radius and the window should be the same.
    IF N_ELEMENTS(RHO) NE N_ELEMENTS(W_FUNC) THEN BEGIN
        MESSAGE, "RHO and W_FUNC have a different number of points"
    ENDIF

    ; Need at least two points to do a numerical integration.
    IF N_ELEMENTS(RHO) LT 2 THEN BEGIN
        MESSAGE, "RHO needs to have at least 2 points"
    ENDIF

    ; Center the independent variable.
    X = (RHO - MEAN(RHO)) / F_SCALE

    ; Step size used for the Riemann sum.
    DRHO = RHO[1] - RHO[0]

    ; The Fresnel kernel, quadratic approximation.
    PSI = 0.5 * !PI * X * X

    ; The complex Fresnel kernel.
    KER = EXP(COMPLEX(0.0, PSI))

    ; Fresnel transform for full transmittance (T = 1).
    T1 = ABS(TOTAL(W_FUNC*KER) * DRHO)

    ; Normalization factor is proportional the reciprocal of the transform.
    NORM_FACT = SQRT(2.0) * F_SCALE / T1
    RETURN, NORM_FACT
END
