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
;       Converts a number to a string.                                         ;
;   Arguments:                                                                 ;
;       NUMBER (real, integer):                                                ;
;           A number.                                                          ;
;   Keywords:                                                                  ;
;       Format (string):                                                       ;
;           Format specifier.                                                  ;
;   Output:                                                                    ;
;       OUT_STRING (string):                                                   ;
;           The input represented as a string.                                 ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Function for converting a number to a string.
FUNCTION NUMBER_TO_STRING, NUMBER, FORMAT = FORMAT

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; Convert the input to a string.
    NUMBER_AS_STRING = STRING(FORMAT = FORMAT, NUMBER)

    ; Remove all unnecessary characters from the resulting string.
    OUT_STRING = STRCOMPRESS(NUMBER_AS_STRING, /REMOVE_ALL)
    RETURN, OUT_STRING
END
