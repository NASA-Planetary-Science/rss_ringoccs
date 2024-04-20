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
;       Computes the second partial derivative of the Fresnel kernel with      ;
;       respect to the azimuth angle Phi.                                      ;
;   Arguments:                                                                 ;
;       R (real, array-like):                                                  ;
;           The ring radius.                                                   ;
;       R0 (real, array-like):                                                 ;
;           The ring radius of the point of interest.                          ;
;       D (real, array-like):                                                  ;
;           Distance from the spacecraft to the ring point.                    ;
;       B (real, array-like):                                                  ;
;           The ring opening angle.                                            ;
;       PHI (real, array-like):                                                ;
;           The ring azimuth angle.                                            ;
;       PHI0 (real, array-like):                                               ;
;           Azimuth angle of the point of interest.                            ;
;   Output:                                                                    ;
;       D2PSI (real, array-like):                                              ;
;           The second derivative d^2 Psi / d Phi^2.                           ;
;   References:                                                                ;
;       1.) Marouf, Tyler, Rosen, 1986                                         ;
;           Profiling Saturn's Rings by Radio Occultation.                     ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;
FUNCTION FRESNEL_KERNEL_D2_PHI, R, R0, D, B, PHI, PHI0

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; Precompute trig values of various parameters.
    COS_B = COS(B)
    COS_PHI0 = COS(PHI0)
    SIN_PHI0 = SIN(PHI0)
    COS_PHI = COS(PHI)
    SIN_PHI = SIN(PHI)

    ; With cos and sin of phi and phi0 computed, the cos and sin of the
    ; difference can be calculated with trig formulas.
    COS_PHI_PHI0 = COS_PHI*COS_PHI0 + SIN_PHI*SIN_PHI0
    SIN_PHI_PHI0 = SIN_PHI*COS_PHI0 - COS_PHI*SIN_PHI0

    ; Save some divisions, precompute the reciprocal of D.
    RCPR_D = 1.0 / D
    RCPR_D_SQ = RCPR_D * RCPR_D

    ; The Xi term from MTR86.
    XI_FACTOR = COS_B * RCPR_D
    XI = XI_FACTOR * (R0*COS_PHI0 - R*COS_PHI)

    ; The Eta term from MTR86.
    ETA = (R0*R0 + R*R - 2.0*R*R0*COS_PHI_PHI0) * RCPR_D_SQ

    ; Several intermediate terms for the second derivative.
    T1 = R * COS_B * COS_PHI * RCPR_D
    T2 = 2.0 * R * R0 * COS_PHI_PHI0 * RCPR_D_SQ
    T3 = R * COS_B * SIN_PHI * RCPR_D
    T4 = 2.0 * R * R0 * SIN_PHI_PHI0 * RCPR_D_SQ
    PSI_FACTOR = 1.0 / SQRT(1.0 + 2.0 * xi + eta)
    PSI_FACTOR_CB = PSI_FACTOR*PSI_FACTOR*PSI_FACTOR
    SQRT_FACTOR = 2.0*T3 + T4
    FACTOR = SQRT_FACTOR*SQRT_FACTOR
    D2PSI_A = (T1 + 0.5*T2) * PSI_FACTOR
    D2PSI_B = 0.25 * FACTOR * PSI_FACTOR_CB
    D2PSI = D2PSI_A - T1 - D2PSI_B
    RETURN, D2PSI
END
