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
;       Calculates the factor in the first iterate of Newton's method          ;
;       for the Fresnel kernel. This is an entirely geometric quantity.        ;
;   Arguments:                                                                 ;
;       R (real, array-like):                                                  ;
;           The radius, dummy variable.                                        ;
;       R0 (real, array-like):                                                 ;
;           Radius of the ring intercept point.                                ;
;       PHI0 (real, array-like):                                               ;
;           The ring azimuth angle.                                            ;
;       B (real, array-like):                                                  ;
;           The ring opening angle.                                            ;
;   Output:                                                                    ;
;       FACTOR (real, array-like):                                             ;
;           The first iterator of Newton's method for the Fresnel scale.       ;
;   References:                                                                ;
;       1.) Marouf, Tyler, Rosen, 1986                                         ;
;           Profiling Saturn's Rings by Radio Occultation.                     ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018/04/10                                                         ;
;------------------------------------------------------------------------------;

; Computes the Newton iterate factor for the Fresnel kernel.
FUNCTION FRESNEL_KERNEL_NEWTON_FACTOR, R, R0, PHI0, B

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Output is scaled by the "normalized" radius.
    SCALE = (R - R0) / R0

    ; Precompute trigonmetric quantities to save on repetition.
    COS_B = COS(B)
    COS_B_SQ = COS_B * COS_b
    COS_P = COS(PHI0)
    SIN_P = SIN(PHI0)
    SIN_P_SQ = SIN_P * SIN_P

    ; The numerator and denominator of the Newton factor.
    NUM = COS_B_SQ * COS_P * SIN_P
    DEN = 1.0 - COS_B_SQ * SIN_P_SQ

    ; The factor is the product of the normalized radius and the ratio of the
    ; two geometric quantities we've computed above. Return this.
    FACTOR = SCALE * NUM / DEN
    RETURN, FACTOR
END
