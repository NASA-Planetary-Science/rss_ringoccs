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
;       Computes the forward Fresnel transform of complex-valued data.         ;
;   Arguments:                                                                 ;
;       T (complex, array-like):                                               ;
;           The complex transmittance.                                         ;
;       KER (complex, array-like):                                             ;
;           The complex-valued Fresnel kernel.                                 ;
;       DX (real):                                                             ;
;           The distance between samples.                                      ;
;       F_SCALE (real):                                                        ;
;           The Fresnel scale.                                                 ;
;   Output:                                                                    ;
;       T_HAT (complex):                                                       ;
;           The Fresnel transform of the input T.                              ;
;   References:                                                                ;
;       1.) Marouf, Tyler, Rosen, 1986                                         ;
;           Profiling Saturn's Rings by Radio Occultation.                     ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018/04/12                                                         ;
;------------------------------------------------------------------------------;
FUNCTION FRESNEL_FORWARD_TRANSFORM, T, KER, DX, F_SCALE

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; The input data should be real or complex valued.
    IF (CHECK_REAL(T) EQ 0) AND (CHECK_COMPLEX(T) EQ 0) THEN BEGIN
        MESSAGE, "T must be real or complex valued."
    ENDIF

    ; The Fresnel kernel should also be real or complex valued.
    IF (CHECK_REAL(KER) EQ 0) AND (CHECK_COMPLEX(KER) EQ 0) THEN BEGIN
        MESSAGE, "KER must be real or complex valued."
    ENDIF

    ; The step size between samples should be a positive real number.
    IF (CHECK_POSITIVE_REAL_NUMBER(DX) EQ 0) THEN BEGIN
        MESSAGE, "DX must be a positive real number."
    ENDIF

    ; Same check for the Fresnel scale.
    IF (CHECK_POSITIVE_REAL_NUMBER(F_SCALE) EQ 0) THEN BEGIN
        MESSAGE, "F_SCALE must be a positive real number."
    ENDIF

    ; Finally, the kernel and the data must be the same size.
    IF N_ELEMENTS(T) NE N_ELEMENTS(KER) THEN BEGIN
        MESSAGE, "T and KER have a different number of elements."
    ENDIF

    ; Scale factor for the integral.
    FACTOR = COMPLEX(0.5, -0.5) * DX / F_SCALE

    ; Approximate using a simple Riemann sum for the integral.
    T_HAT = TOTAL(KER * T) * FACTOR
    RETURN, T_HAT
END
