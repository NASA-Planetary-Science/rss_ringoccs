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
;       Computes the complex-valued Fresnel integrand.                         ;
;   Arguments:                                                                 ;
;       POWER (real, array-like):                                              ;
;           The normalized power.                                              ;
;       FRESNEL_KERNEL (real, array-like):                                     ;
;           The geometric data of the Fresnel kernel.                          ;
;   Output:                                                                    ;
;       T_HAT (complex, array-like):                                           ;
;           The Fresnel integrand from the given geometry and power.           ;
;   References:                                                                ;
;       1.) Marouf, Tyler, Rosen, 1986                                         ;
;           Profiling Saturn's Rings by Radio Occultation.                     ;
;       2.) Goodman, 1969                                                      ;
;           Introduction to Fourier Optics.                                    ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Computes the Fresnel integrand T_HAT of the Fresnel transform.
FUNCTION FRESNEL_INTEGRAND, POWER, FRESNEL_KERNEL

    ; The integrand is ||P|| * e^{-i psi}.
    EXP_MINUS_IPSI = EXP(COMPLEX(0.0, -FRESNEL_KERNEL))
    T_HAT = SQRT(ABS(POWER)) * EXP_MINUS_IPSI
    RETURN, T_HAT
END
