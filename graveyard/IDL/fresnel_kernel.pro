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
;       Computes the real-valued Fresnel kernel.                               ;
;   Arguments:                                                                 ;
;       R (real, array-like):                                                  ;
;           The ring radius.                                                   ;
;       R0 (real, array-like):                                                 ;
;           Dummy variable, ring radius.                                       ;
;       PHI (real, array-like):                                                ;
;           The ring azimuth angle.                                            ;
;       PHI0 (real, array-like):                                               ;
;           Dummy variable, ring azimuth angle.                                ;
;       D (real, array-like):                                                  ;
;           The spacecraft to ring-intercept-point distance.                   ;
;       B (real, array-like):                                                  ;
;           Ring opening angle.                                                ;
;   Output:                                                                    ;
;       PSI (real, array-like):                                                ;
;           The Fresnel kernel from the given geometry.                        ;
;   Notes:                                                                     ;
;       This is the unnormalized Fresnel kernel. It does not include any       ;
;       information about the wavenumber of the incident wave.                 ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Computes the Fresnel kernel from geometric data.
FUNCTION FRESNEL_KERNEL, R, R0, PHI, PHI0, D, B

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; Pre-compute divisions to save time.
    RCPR_D = 1.0 / D
    RCPR_D_SQ = RCPR_D * RCPR_D

    ; The xi and eta factors from the MTR paper.
    XI = (COS(B) * RCPR_D) * (R0*COS(PHI0) - R*COS(PHI))
    ETA = (R0*R0 + R*R - 2.0*R*R0*COS(PHI-PHI0)) * RCPR_D_SQ

    ; The Fresnel kernel as a function of xi and eta.
    PSI = SQRT(1.0 + 2.0*XI + ETA) - (1.0 + XI)
    RETURN, PSI
END
