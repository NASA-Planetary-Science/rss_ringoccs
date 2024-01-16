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
;       Returns the current date as a string in yyyy_mm_dd format.             ;
;   Arguments:                                                                 ;
;       None.                                                                  ;
;   Output:                                                                    ;
;       DATE (string):                                                         ;
;           The current date in yyyy_mm_dd format.                             ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018/04/12                                                         ;
;------------------------------------------------------------------------------;

; Returns the current date in yyyy_mm_dd format.
FUNCTION DATE_STRING

    ; The Julian time can be computed from the built-in SYSTIME function.
    JULIAN_TIME = SYSTIME(/JULIAN)

    ; CALDAT converts Julian time to calendar date.
    CALDAT, JULIAN_TIME, MONTH, DAY, YEAR

    ; Convert the day, month, and year to strings.
    DAY = STRCOMPRESS(STRING(DAY), /REMOVE_ALL)
    MONTH = STRCOMPRESS(STRING(MONTH) ,/REMOVE_ALL)
    YEAR = STRCOMPRESS(STRING(YEAR), /REMOVE_ALL)

    ; For January to September we'll have a 1-character number (0-9) for the
    ; month. To ensure the format is yyyy_mm_dd, prepend a 0.
    IF (STRLEN(month) EQ 1L) THEN BEGIN
        MONTH = '0' + MONTH
    ENDIF

    ; Similarly for the day of the month.
    IF (STRLEN(DAY) EQ 1L) THEN BEGIN
        DAY = '0' + DAY
    ENDIF

    ; We can now create the date string by concatenating.
    DATE = YEAR + '_' + MONTH + '_' + DAY
    RETURN, DATE
END
