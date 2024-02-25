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
;       Computes the normalized equivalent width from a given string.          ;
;   Arguments:                                                                 ;
;       WTYPE (string):                                                        ;
;           The name of the window function.                                   ;
;   Output:                                                                    ;
;       NORM_EQ (real):                                                        ;
;           The normalized equivalent width of the window.                     ;
;   References:                                                                ;
;       1.) Marouf, Tyler, Rosen, 1986                                         ;
;           Profiling Saturn's Rings by Radio Occultation.                     ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018/04/18                                                         ;
;------------------------------------------------------------------------------;

; Returns the normalized equivalent width for pre-computed window functions.
FUNCTION NORMALIZED_EQUIVALENT_WIDTH_FROM_STRING, WTYPE

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; The input should be a string.
    IF IDL_TYPE(WTYPE) NE IDL_TYPE("Hello") THEN BEGIN
        MESSAGE, 'WTYPE must be a string.'
    ENDIF

    ; The normalized equivalent width has been pre-computed for several
    ; different window types. Use these values.
    CASE WTYPE OF
        'rect': NORM_EQ = 1.0
        'coss': NORM_EQ = 1.5
        'kb25': NORM_EQ = 1.65192
        'kb35': NORM_EQ = 1.928446
        'kbmd': NORM_EQ = 1.66

        ; For other windows you should directly compute the value.
        ELSE: MESSAGE, "Invalid window type."
    ENDCASE

    RETURN, NORM_EQ
END
