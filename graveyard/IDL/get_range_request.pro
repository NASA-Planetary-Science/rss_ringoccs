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
;       Creates a range for the rings from a requested range.                  ;
;   Arguments:                                                                 ;
;       RNGREQ (string or array):                                              ;
;           The requested range.                                               ;
;   Output:                                                                    ;
;       RANGE (array):                                                         ;
;           A two element array with the start and end of the range.           ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;
FUNCTION GET_RANGE_REQUEST, RNGREQ

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    INPUT_TYPE = IDL_TYPE(RNGREQ)

    IF INPUT_TYPE NE IDL_TYPE('Hello') AND ~CHECK_REAL(RNGREQ) THEN BEGIN
        MESSAGE, "Input must be a string or a list of real numbers."
    ENDIF

    IF INPUT_TYPE EQ IDL_TYPE('Hello') THEN BEGIN
        REGION = STRCOMPRESS(STRLOWCASE(RNGREQ), /REMOVE_ALL)

        CASE REGION OF
            'maxwell'       : RANGE = [87410.0, 87610.0]
            'maxwellringlet': RANGE = [87410.0, 87610.0]
            'titan'         : RANGE = [77870.0, 77930.0]
            'titanringlet'  : RANGE = [77870.0, 77930.0]
            'huygens'       : RANGE = [117650.0, 117950.0]
            'huygensringlet': RANGE = [117650.0, 117950.0]
            'encke'         : RANGE = [132900.0, 134200.0]
            'enckegap'      : RANGE = [132900.0, 134200.0]
            'all'           : RANGE = [65000.0, 145000.0]
            ELSE: BEGIN
                MESSAGE, "Invalid string for range."
            END
        ENDCASE
    ENDIF ELSE BEGIN

        ; The user should've provided a list or array of numbers.
        IF (N_ELEMENTS(RNGREQ) LT 2) THEN BEGIN
            MESSAGE, "Requested range is a single number."
        ENDIF

        RANGE = DOUBLE(MINMAX(RNGREQ))
    ENDELSE

    RETURN, RANGE
END
