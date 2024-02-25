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
;       Returns a string to a directory and checks for formatting issues.      ;
;   Arguments:                                                                 ;
;       DIRECTORY_STRING (string):                                             ;
;           A directory path.                                                  ;
;   Output:                                                                    ;
;       DIRECTORY (string):                                                    ;
;           The input directory, formatted, and checked for errors.            ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018                                                               ;
;------------------------------------------------------------------------------;

; Returns a string that is the path to a directory and checks for errors.
FUNCTION GET_DIRECTORY, DIRECTORY_STRING

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; Input should be a string, the path to a directory.
    IF IDL_TYPE(DIRECTORY_STRING) NE IDL_TYPE("Hello") THEN BEGIN
        MESSAGE, "Input must be a string"
    ENDIF

    ; Total number of characters in the string.
    LENGTH_OF_STRING = STRLEN(DIRECTORY_STRING)

    ; Ensure the string is properly formatted for whitespace in the name.
    FOR N = 0, LENGTH_OF_STRING - 1 DO BEGIN

        ; Compare characters. If we have a whitespace it should come after an
        ; escape character "\". Check for this.
        CURRENT_CHARACTER = STRMID(DIRECTORY_STRING, N, 1)
        PREVIOUS_CHARACTER = STRMID(DIRECTORY_STRING, N - 1, 1)

        ; Return an error if the string is not properly formatted.
        IF (CURRENT_CHARACTER EQ " ") AND (PREVIOUS_CHARACTER NE "\") THEN BEGIN
            MESSAGE, "'\' must preceed whitespace in directory string."
        ENDIF
    ENDFOR

    ; This is a directory, so it should end with a forward slash.
    LAST_CHARACTER = STRMID(DIRECTORY_STRING, 0, 1, /REVERSE_OFFSET)

    ; If the directory does not end with a forward slash, add it and return.
    IF (LAST_CHARACTER EQ '/') THEN BEGIN
        DIRECTORY = DIRECTORY_STRING
    ENDIF ELSE BEGIN
        DIRECTORY = DIRECTORY_STRING + "/"
    ENDELSE

    RETURN, DIRECTORY
END
