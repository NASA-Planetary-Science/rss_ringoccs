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
;       Computes the Fresnel kernel from precomputed trig values.              ;
;   Arguments:                                                                 ;
;       R (real, array-like):                                                  ;
;           The dummy ring radius.                                             ;
;       R0 (real, array-like):                                                 ;
;           Ring radius of the point.                                          ;
;       D (real, array-like):                                                  ;
;           Distance to the point in the plane.                                ;
;       COSB (real, array-like):                                               ;
;           Cosine of the ring opening angle.                                  ;
;       COSPHI (real, array-like):                                             ;
;           Cosine of the dummy azimuth angle.                                 ;
;       SINPHI (real, array-like):                                             ;
;           Sine of the dummy azimuth angle.                                   ;
;       COSPHI0 (real, array-like):                                            ;
;           Cosine of the azimuth angle.                                       ;
;       SINPHI0 (real, array-like):                                            ;
;           Sine of the azimuth angle.                                         ;
;   Output:                                                                    ;
;       PSI (real, array-like):                                                ;
;           The Fresnel kernel.                                                ;
;   References:                                                                ;
;       1.) Marouf, Tyler, Rosen, 1986                                         ;
;           Profiling Saturn's Rings by Radio Occultation.                     ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Computes the Fresnel kernel from pre-computed trig valued.
FUNCTION FRESNEL_KERNEL_TRIG, R, R0, D, COSB, COSPHI, SINPHI, COSPHI0, SINPHI0

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR,2

    ; Precompute divisions to save a bit of time.
    RCPR_D = 1.0 / D
    RCPR_D_SQ = RCPR_D * RCPR_D

    ; Xi and Eta factors for the kernel, from precomputed trig values.
    XI = COSB * (R0*COSPHI0 - R*COSPHI) * RCPR_D
    ETA = (R0*R0 + R*R - 2.0*R*R0*(SINPHI*SINPHI0+COSPHI*COSPHI0)) * RCPR_D_SQ

    ; Formula for the unnormalized Fresnel kernel.
    PSI = SQRT(1.0 + 2.0 * XI + ETA) - (1.0 + XI)
    RETURN, psi
END
