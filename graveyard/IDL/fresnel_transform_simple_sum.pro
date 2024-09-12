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
;       Simple summation method for the Fresnel transform.                     ;
;   Arguments:                                                                 ;
;       T_IN (array, complex):                                                 ;
;           The complex transmittance.                                         ;
;       KERNEL (array, complex):                                               ;
;           The (complex) Fresnel kernel.                                      ;
;       DX (real):                                                             ;
;           The displacement between points.                                   ;
;       F_SCALE (real):                                                        ;
;           The Fresnel scale. Same units as DX.                               ;
;   Output:                                                                    ;
;       T_OUT (complex):                                                       ;
;           Fresnel transform at the point corresponding to the Fresnel kernel.;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018                                                               ;
;------------------------------------------------------------------------------;
FUNCTION FRESNEL_TRANSFORM_SIMPLE_SUM, T_IN, KERNEL, DX, F_SCALE

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; Fresnel scale must be positive.
    IF ~CHECK_POS_REAL(F_SCALE) THEN BEGIN
        MESSAGE, "Fresnel scale must be a positive real number."
    ENDIF

    ; Similarly, the displacement needs to be a single positive number.
    IF ~CHECK_POS_REAL(DX) THEN BEGIN
        MESSAGE, "Displacement must be a positive real number."
    ENDIF

    ; We naively just do a Riemann sum for the integral. Simple, but
    ; quite often accurate enough.
    RETURN, TOTAL(KERNEL*T_IN) * DX * COMPLEX(1.0,-1.0) / (2.0 * F_SCALE)
END
