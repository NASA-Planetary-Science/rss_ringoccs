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
;       Computes the Fresnel kernel.                                           ;
;   Arguments:                                                                 ;
;       KD (real, array-like):                                                 ;
;           The wavenumber.                                                    ;
;       COSB (real, array-like):                                               ;
;           Cosine of the ring opening angle.                                  ;
;       RHO (real, array-like):                                                ;
;           Ring radius.                                                       ;
;       RHO0 (real, array-like):                                               ;
;           Ring radius for the point of interest.                             ;
;       D (real, array-like):                                                  ;
;           Distance from the ring point to the spacecraft.                    ;
;       PHI0 (real, array-like):                                               ;
;           Ring azimuth angle.                                                ;
;       PHIS (real, array-like):                                               ;
;           Stationary azimuth angle.                                          ;
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
FUNCTION FRESNEL_PSI, KD, COSB, RHO, RHO0, D, PHI0, PHIS

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; Precompute divisions.
    RCPR_D = 1.0 / D
    RCPR_D_SQ = RCPR_D * RCPR_D

    ; Xi and Eta terms from MTR86.
    XI = COSB * (RHO0*COS(PHI0) - RHO*COS(PHIS)) * RCPR_D
    ETA = (RHO0*RHO0 + RHO*RHO - 2.0*RHO*RHO0*COS(PHIS - PHI0)) * RCPR_D_SQ

    ; The Fresnel kernel.
    PSI = KD*(SQRT(1.0 + 2.0*XI + ETA) - (1.0 + XI))
    RETURN, PSI
END
