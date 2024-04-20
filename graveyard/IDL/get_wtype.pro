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
;       Gets the string representing a specific window type.                   ;
;   Arguments:                                                                 ;
;       WINSTR (string):                                                       ;
;           The string for the window type.                                    ;
;   Output:                                                                    ;
;       WTYPE (string):                                                        ;
;           The window type, in lowercase and with white space removed.        ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018/04/04                                                         ;
;------------------------------------------------------------------------------;
FUNCTION GET_WTYPE, WINSTR

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; Input must be a string. Check for this.
    IF IDL_TYPE(WINSTR) NE IDL_TYPE("Hello") THEN BEGIN
        MESSAGE, "Input must be a string."
    ENDIF

    ; Legal strings for window types.
    TRECT = ["rect", "rec", "rectangle", "rectangular"]
    TCOSS = ["coss", "cos", "cosine", "raisedcosine", "cossq", "sqcos"]
    TKB25 = ["kb25", "k25", "kb_25", "kb2.5", "kb_2.5", "kaiserbessel_25"]
    TKB35 = ["kb35", "k35", "kb_35", "kb3.5", "kb_3.5", "kaiserbessel_35"]
    TKBMD = ["kbmd", "kmd", "kbmod", "kb_mod", "kaiserbessel_mod", "kbmodified"]

    ; Remove white space and make lower case.
    WSTR  = STRCOMPRESS(STRLOWCASE(WINSTR), /REMOVE_ALL)

    ; Check if the input string matches a valid window type.
    IF WHERE(TRECT EQ WSTR) GE 0 THEN BEGIN
        RETURN, "rect"
    ENDIF

    IF WHERE(TCOSS EQ WSTR) GE 0 THEN BEGIN
        RETURN, "coss"
    ENDIF

    IF WHERE(TKB25 EQ WSTR) GE 0 THEN BEGIN
        RETURN, "kb25"
    ENDIF

    IF WHERE(TKB35 EQ WSTR) GE 0 THEN BEGIN
        RETURN, "kb35"
    ENDIF

    IF WHERE(TKBMD EQ WSTR) GE 0 THEN BEGIN
        RETURN, "kbmd"
    ENDIF

    MESSAGE, 'Illegal WTYPE. Allowed: "rect","coss","kb25","kb35","kbmd"'
END
