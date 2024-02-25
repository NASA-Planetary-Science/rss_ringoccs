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
;       Experiments with Fresnel diffraction for a square well.                ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Tells the compiler that integers should be 32 bits, not 16.
COMPILE_OPT IDL2

; Error checking code.
ON_ERROR, 2

; User inputs for the experiment.
LENGTH = 20.0   ; Length of data set, in km
N = 6E4         ; Number of data points
M = 5E3         ; Number of points in W plot
HEIGHT = 1.0    ; Height of well
DEPTH = 0.5     ; Depth of well
WIDTH = 1.0     ; Width of well
LAMBDA = 1.0E-6 ; Wavelength of light
D = 1.0E5       ; Distance from observer to well.
P1 = 0.7        ; Parameter for W plot
P2 = 0.6        ; Parameter for W plot

; Computed values from those inputs.
HALF_N = N / 2
QUARTER_N = N / 4
X = FINDGEN(N) / FLOAT(CEIL((N / LENGTH))) - 0.5*LENGTH
ABS_X = ABS(X)
F_RAN = dblarr(M)
DX = x[1] - x[0]
F = sqrt(LAMBDA * D * 0.5)
RCPR_F = 1.0 / F
A = -0.5*WIDTH
B = 0.5*WIDTH
FRESNEL_FACTOR = -0.5 * !PI * (X * X) / (F * F)
ARG = EXP(COMPLEX(0.0, FRESNEL_FACTOR))
I_ARG = EXP(COMPLEX(0.0, -FRESNEL_FACTOR))
T_HAT_IDEAL = FRESNEL_TRANSFORM_SQUARE_WELL(X, A, B, F, HEIGHT, DEPTH)
FFT_T_HAT_IDEAL = SHIFT(FFT(T_HAT_IDEAL, 1), HALF_N)
P_IDEAL = ABS(T_HAT_IDEAL)
CHARACTER = FLTARR(M)
WINDOW_RANGE_FACTOR = ROUND(WIDTH / DX)
LEFT_WINDOW_RANGE = N/2 - WINDOW_RANGE_FACTOR + 1
RIGHT_WINDOW_RANGE = N/2 + WINDOW_RANGE_FACTOR - 1
FORWARD_FACTOR = COMPLEX(0.5, 0.5) * RCPR_F
INVERSE_FACTOR = COMPLEX(0.5, -0.5) * RCPR_F
T = FLTARR(N) + HEIGHT
T[LEFT_WINDOW_RANGE:RIGHT_WINDOW_RANGE] = DEPTH
FFT_T = SHIFT(FFT(T, 1), HALF_N)

; Loop through different window sizes.
FOR N_ITER = 1, M DO BEGIN
    WINDOW_PERCENT = 1.0 * (FLOAT(N_ITER) / FLOAT(M))
    WINDOW_WIDTH = WINDOW_PERCENT * LENGTH
    X_RANGE = WHERE(ABS_X LE 0.5 * WINDOW_WIDTH)
    W = FLTARR(N)
    W[X_RANGE] = TRIANGULAR_WINDOW(X[X_RANGE])
    W_ARG = ARG * W
    WI_ARG = I_ARG * W
    FFT_ARG = SHIFT(FFT(W_ARG, 1), HALF_N)
    FFT_I_ARG = SHIFT(FFT(WI_ARG, 1), HALF_N)
    KERN = FFT_T * FFT_ARG
    INV_KERN = FFT_T_HAT_IDEAL * FFT_I_ARG
    T_HAT = INVERSE_FACTOR * SHIFT(FFT(KERN, -1), HALF_N)
    P = DX * ABS(T_HAT)
    T_IDEAL_RECON = FORWARD_FACTOR * SHIFT(FFT(INV_KERN, -1), HALF_N)
    P_IDEAL_RECON = DX * ABS(T_IDEAL_RECON)
    P_IDEAL_BASE = MEAN(P_IDEAL_RECON[100:QUARTER_N])
    ALPHA1 = MIN(WHERE(P_IDEAL_RECON LT P1*P_IDEAL_BASE))
    ALPHA2 = MIN(WHERE(P_IDEAL_RECON LT P2*P_IDEAL_BASE))
    CHAR = FLOAT((ALPHA2 - ALPHA1)) * RCPR_F
    CHARACTER[N_ITER - 1] = CHAR
    F_RAN[N_ITER - 1] = P_IDEAL_BASE
    PLOT, X, P_IDEAL_RECON, YRANGE = [-0.1, 1.4]
    OPLOT, X, W, COLOR = CGCOLOR('red')
    WAIT, 1.0E-3
    PRINT, N_ITER, P_IDEAL_BASE
ENDFOR

PLOT, X, P_IDEAL
OPLOT, X, P, COLOR = CGCOLOR('red')
OPLOT, X, T, COLOR = CGCOLOR('green')
OPLOT, X, P_IDEAL_RECON, COLOR = CGCOLOR('yellow')
OPLOT, X, 0.5 * W, COLOR = CGCOLOR('blue')

END
