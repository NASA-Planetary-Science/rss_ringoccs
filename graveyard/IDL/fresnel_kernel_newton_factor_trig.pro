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
;       Computes the Newton iterate for the Fresnel kernel from pre-computed   ;
;       trigonometric values.                                                  ;
;   Arguments:                                                                 ;
;       R (real, array-like):                                                  ;
;           The dummy ring radius.                                             ;
;       R0 (real, array-like):                                                 ;
;           Ring radius of the point.                                          ;
;       COSB (real, array-like):                                               ;
;           Cosine of the ring opening angle.                                  ;
;       COSPHI0 (real, array-like):                                            ;
;           Cosine of the azimuth angle.                                       ;
;       SINPHI0 (real, array-like):                                            ;
;           Sine of the azimuth angle.                                         ;
;   Output:                                                                    ;
;       FACTOR (real, array-like):                                             ;
;           The Newton iterate for the Fresnel kernel.                         ;
;   References:                                                                ;
;       1.) Marouf, Tyler, Rosen, 1986                                         ;
;           Profiling Saturn's Rings by Radio Occultation.                     ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Computes the Newton iterate factor for the Fresnel kernel.
FUNCTION FRESNEL_KERNEL_NEWTON_FACTOR_TRIG, R, R0, COSB, COSPHI0, SINPHI0

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; Precompute trigonmetric quantities to save on repetition.
    COSB_SQ = COSB * COSB
    SINPHI0_SQ = SINPHI0 * SINPHI0

    ; Numerator and denominator for the Newton iterate.
    NUM = (R - R0) * COSB_SQ*COSPHI0*SINPHI0
    DEN = (1.0 - COSB_SQ*SINPHI0_SQ) * R0

    ; Compute and return.
    FACTOR = NUM / DEN
    RETURN, FACTOR
END
