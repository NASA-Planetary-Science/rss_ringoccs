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
PRO WINDOW_CALCULATION, PARAMETERS, IN_STRUCT, W_KM_VALS

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; Compute the normalized equivalent width with respect to window type.
    CASE PARAMETERS.WTYPE OF
        'rect': NORM_EQ = 1.0
        'cos' : NORM_EQ = 1.50
        'kb25': NORM_EQ = 1.6519200
        'kb35': NORM_EQ = 1.9284460
        'kbmod': NORM_EQ = 1.660
    ENDCASE

    ; From USO presentation.
    SIGMA = 2.0E-13

    ; Angular frequency (Hz).
    OMEGA = 2.0 * !PI * IN_STRUCT.F_SKY_HZ_VALS

    ; Parameters for inverse of Resolution = f(Window).
    FACTOR = OMEGA * SIGMA
    ALPHA = 0.5 * FACTOR * FACTOR / IN_STRUCT.RHO_DOT_KMS_VALS
    P = PARAMETERS.RES_KM / (ALPHA * IN_STRUCT.F_KM_VALS * IN_STRUCT.F_KM_VALS)

    ; The inverse exists only if P > 1.
    IF MIN(P) LE 1.0001 THEN BEGIN
        MESSAGE, 'ERROR: Either rho_dot, the f_scale, or res_km is too small.'
    ENDIF

    ; Create window variable, window width (in km) for each point.
    W_KM_VALS = SQRT(2.0) * NORM_EQ * INVERSE_RESOLUTION(P) / ALPHA

    ; If the window width is large, inversion will take a while.
    IF MAX(W_KM_VALS GT 1000.0) THEN BEGIN
        MESSAGE, 'Max window width required is too large.'
    ENDIF
END
