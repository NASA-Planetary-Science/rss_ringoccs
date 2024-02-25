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
;       Determines the type number of an IDL data type.                        ;
;   Arguments:                                                                 ;
;       VAR (Any type):                                                        ;
;           Any input.                                                         ;
;   Output:                                                                    ;
;       TYPE_NUMBER (Int):                                                     ;
;           Integer that represents the IDL data type of the input.            ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018/04/18                                                         ;
;------------------------------------------------------------------------------;

; Function for determining the data type of an input.
FUNCTION IDL_TYPE, VAR

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; The size function has the data type in the second-to-last entry.
    SIZE_VAR = SIZE(VAR)

    ; IDL supports negative indexing. Return the second-to-last entry.
    RETURN, SIZE_VAR[-2]
END
