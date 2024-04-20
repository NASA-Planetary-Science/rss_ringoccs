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
;       Computes the 2nd derivative of the Fresnel kernel with respect to phi. ;
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
;       D2PSI (real, array-like):                                              ;
;           The derivative d^2 Psi / d Phi^2.                                  ;
;   References:                                                                ;
;       1.) Marouf, Tyler, Rosen, 1986                                         ;
;           Profiling Saturn's Rings by Radio Occultation.                     ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;
FUNCTION FRESNEL_PSI_XX, KD, COSB, RHO, RHO0, D, PHI0, PHIS

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; Precompute divisions.
    RCPR_D = 1.0 / D
    RCPR_D_SQ = RCPR_D * RCPR_D

    ; Precompute trig values.
    SIN_PHI0 = SIN(PHI0)
    COS_PHI0 = COS(PHI0)
    SIN_PHIS = SIN(PHIS)
    COS_PHIS = COS(PHIS)

    ; With cos and sin of phi0 and phi_s computed, the cos and sin of the
    ; difference can be calculated with trig formulas.
    COS_PHIS_PHI0 = COS_PHIS*COS_PHI0 + SIN_PHIS*SIN_PHI0
    SIN_PHIS_PHI0 = SIN_PHIS*COS_PHI0 - COS_PHIS*SIN_PHI0

    ; Xi and Eta terms from MTR86.
    XI = COSB * (RHO0 * COS_PHI0 - RHO * COS_PHIS) * RCPR_D
    ETA = (RHO*RHO + RHO0*RHO0 - 2.0*RHO*RHO0*COS_PHIS_PHI0) * RCPR_D_SQ

    ; Derivative factors for Xi and Eta.
    XI_X = RHO * COSB * SIN_PHIS * RCPR_D
    ETA_X = 2.0 * RHO * RHO0 * SIN_PHIS_PHI0 * RCPR_D_SQ

    XI_XX = RHO * COSB * COS_PHIS * RCPR_D
    ETA_XX = 2.0 * RHO * RHO0 * COS_PHIS_PHI0 * RCPR_D_SQ

    ; Fresnel kernel factor.
    PSI0 = SQRT(1.0 + 2.0*XI + ETA)
    FRESNEL_FACTOR = 1.0 / PSI0
    FRESNEL_FACTOR_CB = FRESNEL_FACTOR*FRESNEL_FACTOR*FRESNEL_FACTOR
    SQRT_FACTOR = 2.0*ETA_X + XI_X
    FACTOR = SQRT_FACTOR*SQRT_FACTOR

    D2PSI_A = (eta_xx + 0.5*xi_xx) * FRESNEL_FACTOR
    D2PSI_B = 0.25 * FACTOR * FRESNEL_FACTOR_CB

    ; The second derivative d^2 Psi / d Phi^2.
    D2PSI = D2PSI_A - D2PSI_B - XI_XX
    RETURN, KD * D2PSI
END
