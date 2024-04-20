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
;       COSB (real, array-like):                                               ;
;           Cosine of the ring opening angle.                                  ;
;       COSPHI (real, array-like):                                             ;
;           Cosine of the ring azimuth angle.                                  ;
;       SINPHI (real, array-like):                                             ;
;           Sine of the ring azimuth angle.                                    ;
;       COSPHI0 (real, array-like):                                            ;
;           Cosine of the azimuth angle of the point of interest.              ;
;       SINPHI0 (real, array-like):                                            ;
;           Sine of the azimuth angle of the point of interest.                ;
;   Output:                                                                    ;
;       D2PSI (real, array-like):                                              ;
;           The second derivative d^2 Psi / d Phi^2.                           ;
;   References:                                                                ;
;       1.) Marouf, Tyler, Rosen, 1986                                         ;
;           Profiling Saturn's Rings by Radio Occultation.                     ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018/04/10                                                         ;
;------------------------------------------------------------------------------;
FUNCTION FRESNEL_KERNEL_D2_PHI_TRIG, R, R0, D, COSB, COSPHI, SINPHI, COSPHI0, SINPHI0

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; Save on divisions by precomputed reciprocals.
    RCPR_D = 1.0 / D
    RCPR_D_SQ = RCPR_D * RCPR_D

    ; cos(x - y) and sin(x - y) can be computed from the precomputed values
    ; cos(x), cos(y), sin(y), and sin(y). No need for additional sin and cos.
    COS_PHI_PHI0 = SINPHI*SINPHI0 + COSPHI*COSPHI0
    SIN_PHI_PHI0 = SINPHI*COSPHI0 - SINPHI0*COSPHI

    ; Xi term from MTR86.
    XI_FACTOR = COSB * RCPR_D
    XI = XI_FACTOR * (R0*COSPHI0 - R*COSPHI)

    ; Eta term from MTR86.
    ETA = (R0*R0 + R*R - 2.0*R*R0*COS_PHI_PHI0) * RCPR_D_SQ

    ; Lots of intermediate terms for the first and second derivatives.
    T1 = R * COSB * COSPHI * RCPR_D
    T2 = 2.0 * R * R0 * COS_PHI_PHI0 * RCPR_D_SQ
    T3 = R * COSB * SINPHI * RCPR_D
    T4 = 2.0 * R * R0 * SIN_PHI_PHI0 * RCPR_D_SQ
    PSI_FACTOR = 1.0 / SQRT(1.0 + 2.0 * XI + ETA)
    PSI_FACTOR_CB = PSI_FACTOR * PSI_FACTOR * PSI_FACTOR
    SQRT_FACTOR = 2.0*T3 + T4
    FACTOR = SQRT_FACTOR * SQRT_FACTOR
    D2PSI_A = (T1 + 0.5*T2)*PSI_FACTOR
    D2PSI_B = 0.25 * FACTOR * PSI_FACTOR_CB
    D2PSI = D2PSI_A - T1 - D2PSI_B
    RETURN, D2PSI
END
