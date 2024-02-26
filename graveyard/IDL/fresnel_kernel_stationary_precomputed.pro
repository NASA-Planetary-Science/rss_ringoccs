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
;       Computes the stationary Fresnel kernel using precomputed values.       ;
;   Arguments:                                                                 ;
;       KD (real, array-like):                                                 ;
;           The product of the wave number and distance to ring plane.         ;
;       COS_B (real, array-like):                                              ;
;           The cosine of the ring opening angle B.                            ;
;       R (real, array-like):                                                  ;
;           Dummy variable, ring radius.                                       ;
;       R0 (real, array-like):                                                 ;
;           Radius of the point of interest.                                   ;
;       D (real, array-like):                                                  ;
;           The spacecraft to ring-intercept-point distance.                   ;
;       PHI0 (real, array-like):                                               ;
;           The ring azimuth angle.                                            ;
;       PHI_S (real, array-like):                                              ;
;           The pre-computed stationary azimuth angle.                         ;
;   Output:                                                                    ;
;       PSI (real, array-like):                                                ;
;           The Fresnel kernel from the given geometry.                        ;
;   References:                                                                ;
;       1.) Marouf, Tyler, Rosen, 1986                                         ;
;           Profiling Saturn's Rings by Radio Occultation.                     ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Computes the stationary Fresnel kernel from pre-computed stationary angles.
FUNCTION FRESNEL_KERNEL_STATIONARY_PRECOMPUTED, KD, COS_B, R, R0, D, PHI0, PHI_S

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; Array for the output Fresnel kernel.
    PSI_VALS = FLTARR(5, N_ELEMENTS(R))

    ; Precompute a bunch of values to save on time.
    COS_PHI0 = COS(PHI0)
    RCPR_D = 1.0 / D
    RCPR_D_SQ = RCPR_D * RCPR_D
    R_SQ = R*R
    R0_SQ = R0*R0
    FACTOR = 2.0 * R0 * R

    ; Compute the Fresnel kernel for each given stationary value.
    FOR N = 0, 4 DO BEGIN
        XI = COS_B * (R0 * COS_PHI0 - R * COS(PHI_S[N])) * RCPR_D
        ETA = (R0_SQ + R_SQ - FACTOR * COS(PHI_S[N] - PHI0)) * RCPR_D_SQ
        PSI_VALS[N, *] = KD * (SQRT(1.0 + 2.0*XI + ETA) - (1.0 + XI))
    ENDFOR

    RETURN, PSI_VALS
END
