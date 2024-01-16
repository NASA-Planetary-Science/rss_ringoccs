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
;       Computes the normalized sinc function.                                 ;
;   Arguments:                                                                 ;
;       VAR (real, complex, array-like):                                       ;
;           The independent variable.                                          ;
;   Output:                                                                    ;
;       SINC_X (real, complex, array-like):                                    ;
;           The normalized sinc of the input.                                  ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Function for computing the normalized sinc function.
FUNCTION SINC, VAR

    ; Error checking code.
    ON_ERROR, 2

    ; We allow real and complex inputs.
    IS_REAL = CHECK_REAL(VAR)
    IS_COMPLEX = CHECK_COMPLEX(VAR)

    ; Number of elements in the input array.
    VAR_LENGTH = N_ELEMENTS(VAR)

    ; If the input is neither real nor complex, abort.
    IF (IS_REAL EQ 0) AND (IS_COMPLEX EQ 0) THEN BEGIN
        MESSAGE, "Input should be real or complex."
    ENDIF

    ; There should be at least one element in the input array.
    IF (VAR_LENGTH LT 1) THEN BEGIN
        MESSAGE, "Input is undefined."
    ENDIF

    ; Get the values where the array is non-zero.
    NON_ZERO = WHERE(VAR NE 0.0)

    ; Create an array and set all entries to 1. Preserve the type of the input.
    IF (IS_REAL EQ 1) THEN BEGIN
        SINC_VAR = REPLICATE(1.0, VAR_LENGTH)
    ENDIF ELSE BEGIN
        SINC_VAR = REPLICATE(COMPLEX(1.0, 0.0), VAR_LENGTH)
    ENDELSE

    ; For all non-zero values use the formula sinc(z) = sin(pi z) / (pi z).
    IF (MIN(NON_ZERO) NE -1) THEN BEGIN
        SINC_VAR[NON_ZERO] = SIN(!PI * VAR[NON_ZERO]) / (!PI * VAR[NON_ZERO])
    ENDIF

    RETURN, SINC_VAR
END
