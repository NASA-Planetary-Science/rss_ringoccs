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
;       Determines if a data type is complex-valued.                           ;
;   Arguments:                                                                 ;
;       VAR (Any type):                                                        ;
;           Any input.                                                         ;
;   Output:                                                                    ;
;       IS_COMPLEX (Int):                                                      ;
;           Integer that determines if the input is complex-valued.            ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018/04/18                                                         ;
;------------------------------------------------------------------------------;

; Function for checking if a variable is complex-valued.
FUNCTION CHECK_COMPLEX, VAR

    ; The type of the input.
    VAR_TYPE = IDL_TYPE(VAR)

    ; The type numbers for complex variables.
    COMPLEX_TYPES = [IDL_TYPE(COMPLEX(0)), IDL_TYPE(DCOMPLEX(0))]

    ; Check if we have a match. If not, this is not complex-valued.
    IF (WHERE(VAR_TYPE EQ COMPLEX_TYPES) EQ -1) THEN BEGIN
        RETURN, 0
    ENDIF

    ; Otherwise we have a complex variable. Return true.
    RETURN, 1
END
