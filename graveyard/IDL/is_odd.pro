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
;       Determines if a variable is an odd integer.                            ;
;   Arguments:                                                                 ;
;       VAR (Any type):                                                        ;
;           Any input.                                                         ;
;   Output:                                                                    ;
;       ODD (Int):                                                             ;
;           An integer, 0 or 1, depending on whether the input is odd.         ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Function for determining if an input is an odd number.
FUNCTION IS_ODD, VAR

    ; The data type of the input. Should be an integer type.
    VAR_TYPE = IDL_TYPE(VAR)

    ; Total number of elements in the input. Should be one for a number.
    VAR_LENGTH = N_ELEMENTS(VAR)

    ; The type numbers for the IDL integer data types.
    INTEGER_TYPES = [               $
        IDL_TYPE(0),                $
        IDL_TYPE(LONG(0)),          $
        IDL_TYPE(UINT(0)),          $
        IDL_TYPE(ULONG(0)),         $
        IDL_TYPE(LONG64(0)),        $
        IDL_TYPE(ULONG64(0))        $
    ]

    ; Odd numbers should not be array objects, just a single integer.
    IF VAR_LENGTH NE 1L THEN BEGIN
        RETURN, 0
    ENDIF

    ; Odd numbers are represented by integer data types.
    IF (WHERE(VAR_TYPE EQ INTEGER_TYPES) EQ -1) THEN BEGIN
        RETURN, 0
    ENDIF

    ; Parity is the same for positive and negative numbers.
    ABS_VAR = ABS(VAR)

    ; Odd numbers are multiples of two plus one.
    IF (2L * (ABS_VAR / 2L) EQ ABS_VAR) THEN BEGIN
        RETURN, 0
    ENDIF

    ; All checks passed, we have an odd number.
    RETURN, 1
END
