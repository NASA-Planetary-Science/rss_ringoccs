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
;       Determines if a data type is a positive real number.                   ;
;   Arguments:                                                                 ;
;       VAR (Any type):                                                        ;
;           Any input.                                                         ;
;   Output:                                                                    ;
;       IS_POS_REAL (Int):                                                     ;
;           Integer that determines if the input is a positive real number.    ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018/04/12                                                         ;
;------------------------------------------------------------------------------;

; Function for determining if this input is a positive real number.
FUNCTION CHECK_POSITIVE_REAL_NUMBER, VAR

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; The type of the input.
    VAR_TYPE = IDL_TYPE(VAR)

    ; The number of elements the variable represents.
    VAR_LENGTH = N_ELEMENTS(VAR)

    ; Data types corresponding to real numbers.
    REAL_TYPES = [                  $
        IDL_TYPE(0),                $
        IDL_TYPE(LONG(0)),          $
        IDL_TYPE(UINT(0)),          $
        IDL_TYPE(ULONG(0)),         $
        IDL_TYPE(LONG64(0)),        $
        IDL_TYPE(ULONG64(0)),       $
        IDL_TYPE(0.0),              $
        IDL_TYPE(0.d0)              $
    ]

    ; If the data type does not match, we do not have something real.
    IF (WHERE(VAR_TYPE EQ REAL_TYPES) EQ -1) THEN BEGIN
        RETURN, 0
    ENDIF

    ; We are looking for a single number. The length should be 1.
    IF VAR_LENGTH NE 1 THEN BEGIN
        RETURN, 0
    ENDIF

    ; If we have a negative or zero-valued input, return false.
    IF VAR LE 0.0 THEN BEGIN
        RETURN, 0
    ENDIF

    ; Otherwise, we have a positive real-valued number. Return true.
    RETURN, 1
END
