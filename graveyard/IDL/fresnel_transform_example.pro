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
;       Creates a plot of the Fresnel transform of a simple square well.       ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Tells the compiler that integers should be 32 bits, not 16.
COMPILE_OPT IDL2

; Number of points for the x-axis.
NUMBER_OF_POINTS = 1024
HALF_NUMBER_OF_POINTS = NUMBER_OF_POINTS / 2

; Factor for the inverse Fourier transform.
FACTOR = SQRT(!PI * 0.5) / FLOAT(NUMBER_OF_POINTS)

; Set spatial resolution.
DX = 0.01

; Independent variable, and the square of it.
U_N = (FINDGEN(NUMBER_OF_POINTS) - NUMBER_OF_POINTS*0.5)*DX
U_N_SQ = U_N * U_N

; Create an impulse array for Fresnel transform.
U = (FLTARR(NUMBER_OF_POINTS) - NUMBER_OF_POINTS * 0.5)*DX

; Create an impulse profile.
U[HALF_NUMBER_OF_POINTS - 5:HALF_NUMBER_OF_POINTS + 5] = 1.0

; Set z value.
Z = 2.0

; Choose wavelength to be used.
LAMBDA = 1.0

; Compute the Frensel Scale.
F_SQ = LAMBDA * Z * 0.5

; Perform Fourier transform.
U_HAT = SHIFT(FFT(U, 1), HALF_NUMBER_OF_POINTS)

; Create the kernel function, e^{i psi} = e^{i pi (x / 2F)^2}.
FRESNEL_FACTOR = !PI / (2.0 * F_SQ)
H = EXP(COMPLEX(0.0, FRESNEL_FACTOR * U_N_SQ))

; Fourier transform of the kernel function.
H_HAT = SHIFT(FFT(H, 1), HALF_NUMBER_OF_POINTS)

; Convolution theorem: Fourier transform of a convolution is the point-wise
; product of the individual Fourier transforms.
OH = U_HAT * H_HAT

; Compute the inverse Fourier transform.
T_HAT = SHIFT(FFT(OH, -1), HALF_NUMBER_OF_POINTS) * FACTOR

; Compute the power
SQRT_POWER = ABS(T_HAT)
POWER = SQRT_POWER * SQRT_POWER

; Plot the results.
PLOT, U_N, POWER
OPLOT, U_N, U

END
