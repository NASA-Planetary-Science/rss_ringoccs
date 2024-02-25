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
;       Computes the diffraction pattern of a square well.                     ;
;   Arguments:                                                                 ;
;       X (real, array-like):                                                  ;
;           The independent variable.                                          ;
;       A (real):                                                              ;
;           Left end point of the square well.                                 ;
;       B (real):                                                              ;
;           Right end point of the square well.                                ;
;       F (real):                                                              ;
;           Fresnel scale.                                                     ;
;       HEIGHT (real):                                                         ;
;           Height of the square well.                                         ;
;       WIDTH (real):                                                          ;
;           Width of the square well.                                          ;
;   Output:                                                                    ;
;       T_HAT (complex, array-like):                                           ;
;           The diffraction pattern of a square well.                          ;
;   References:                                                                ;
;       1.) Goodman, 1969                                                      ;
;           Introduction to Fourier Optics.                                    ;
;       2.) Carroll, Ostie, 2007                                               ;
;           An Introduction to Modern Astrophysics.                            ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Computes the Fresnel transform for a square well.
FUNCTION FRESNEL_TRANSFORM_SQUARE_WELL, X, A, B, F, HEIGHT, DEPTH

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; Constant factors in the formula.
    RCPR_F = 1.0 / F
    LEFT_SHIFT = (A - X) * RCPR_F
    RIGHT_SHIFT = (B - X) * RCPR_F

    ; Real and imaginary parts for the transform.
    REAL_PART = FRESNEL_COS(RIGHT_SHIFT) - FRESNEL_COS(LEFT_SHIFT)
    IMAG_PART = FRESNEL_SIN(RIGHT_SHIFT) - FRESNEL_SIN(LEFT_SHIFT)

    ; Scale factor for the integral.
    SCALE = (HEIGHT - DEPTH) * COMPLEX(0.5, -0.5)

    ; The diffraction pattern.
    T_HAT = HEIGHT - SCALE * COMPLEX(REAL_PART, IMAG_PART)
    RETURN, T_HAT
END
