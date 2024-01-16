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
;       Determines if a data type is a ternary value.                          ;
;   Arguments:                                                                 ;
;       VAR (Any type):                                                        ;
;           Any input.                                                         ;
;   Output:                                                                    ;
;       IS_TERNARY (Int):                                                      ;
;           Integer that determines if the input is a ternary value.           ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018/04/18                                                         ;
;------------------------------------------------------------------------------;
FUNCTION CHECK_TERNARY, VAR

    ; The data type of the input. Should be an integer type.
    VAR_TYPE = IDL_TYPE(VAR)

    ; Total number of elements in the input. Should be one for a ternary value.
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

    ; A ternary should not be an array object, just a single integer.
    IF VAR_LENGTH NE 1L THEN BEGIN
        RETURN, 0
    ENDIF

    ; Ternary values are represented by integer data types.
    IF (WHERE(VAR_TYPE EQ INTEGER_TYPES) EQ -1) THEN BEGIN
        RETURN, 0
    ENDIF

    ; Ternary values are zero, one, or two.
    IF (VAR LT 0) OR (VAR GT 2) THEN BEGIN
        RETURN, 0
    ENDIF

    ; All checks passed, we most likely have a ternary input.
    RETURN, 1
END
