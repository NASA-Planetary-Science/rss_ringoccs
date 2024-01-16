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
;       Computes the complex-valued Fresnel sine.                              ;
;   Arguments:                                                                 ;
;       VAR (real, complex, array-like):                                       ;
;           The independent variable.                                          ;
;   Output:                                                                    ;
;       F_SIN (real, complex, array-like):                                     ;
;           The Fresnel sine of the input.                                     ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Computes the Fresnel sine of the input.
FUNCTION FRESNEL_SIN, VAR

    ; Error checking code.
    ON_ERROR, 2

    ; Check for a valid input. Must have at least one element.
    IF N_ELEMENTS(VAR) EQ 0L THEN BEGIN
        MESSAGE, 'Input is undefined.'
    ENDIF

    ; This function allows real and complex arrays.
    IS_REAL = CHECK_REAL(VAR)
    IS_COMPLEX = CHECK_COMPLEX(VAR)

    ; Any other type of input is invalid. Raise an error.
    IF (IS_REAL EQ 0L) AND (IS_COMPLEX EQ 0L) THEN BEGIN
        MESSAGE, 'Input must be real or complex valued'
    ENDIF

    ; The normalized Fresnel functions have factors of root pi in the
    ; the argument, whereas the Error function does not. We'll need to scale
    ; the input by this value.
    SCALE_FACTOR = VAR * SQRT(!pi) / 2.0

    ; The two direction in the complex plane used to scale the argument.
    NORTH_EAST = COMPLEX(1.0, 1.0)
    SOUTH_EAST = COMPLEX(1.0, -1.0)

    ; The Fresnel sine is a sum of two error functions. Compute these.
    LEFT_TERM = (NORTH_EAST / 4.0) * ERF(NORTH_EAST * SCALE_FACTOR)
    RIGHT_TERM = (SOUTH_EAST / 4.0) * ERF(SOUTH_EAST * SCALE_FACTOR)

    F_SIN = LEFT_TERM + RIGHT_TERM

    ; If the input is real, we should return the output as a real variable.
    IF IS_REAL THEN BEGIN
        RETURN, REAL_PART(F_SIN)
    ENDIF

    ; Otherwise a complex output is suitable.
    RETURN, F_SIN
END
