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
;       Creates the string for a save file of a given struct.                  ;
;   Arguments:                                                                 ;
;       STRUCT (structure):                                                    ;
;           Structure with the reconstruction data.                            ;
;   Keywords:                                                                  ;
;       DIR (string):                                                          ;
;           Path to the save directory.                                        ;
;       REV (string):                                                          ;
;           Rev number for the data.                                           ;
;       LAGRANGE (ternary):                                                    ;
;           Mode of numerical differentiation used.                            ;
;       NORMALIZE (bool):                                                      ;
;           Boolean for normalizing the window width.                          ;
;   Output:                                                                    ;
;       OUT_STRING (string):                                                   ;
;           The file name for the data.                                        ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018/04/12                                                         ;
;------------------------------------------------------------------------------;
FUNCTION GET_OUT_STRING, STRUCT, DIR = DIR, REV = REV, $
                         LAGRANGE = LAGRANGE, NORMALIZE = NORMALIZE

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    IF KEYWORD_SET(DIR) THEN BEGIN

        ; The directory should be a path, and hence a string.
        IF IDL_TYPE(DIR) NE IDL_TYPE("Hello") THEN BEGIN
            MESSAGE, "DIR must be a string."
        ENDIF

        DIR = GET_DIRECTORY(DIR)

    ENDIF ELSE BEGIN
        DIR = ""
    ENDELSE

    IF KEYWORD_SET(REV) THEN BEGIN

        ; Rev should be a string with information about the occultation.
        IF IDL_TYPE(REV) NE IDL_TYPE("Hello") THEN BEGIN
            MESSAGE, "REV must be a string."
        ENDIF

        REV = STRCOMPRESS(REV, /REMOVE_ALL)

    ENDIF ELSE BEGIN
        REV = ""
    ENDELSE

    IF KEYWORD_SET(LAGRANGE) THEN BEGIN

        ; Three options for differentiation, normal, Lagrange, and fast.
        ; This variable should be a ternary value.
        IF ~CHECK_TERNARY(LAGRANGE) THEN BEGIN
            MESSAGE, "LAGRANGE must be Ternary."
        ENDIF

        LAGRANGE = LONG(LAGRANGE)

    ENDIF ELSE BEGIN
        LAGRANGE = 0L
    ENDELSE

    IF KEYWORD_SET(NORMALIZE) THEN BEGIN

        ; Normalize is a Boolean, whether or not to normalize the output
        ; by the window width.
        IF ~CHECK_BOOLE(NORMALIZE) THEN BEGIN
            MESSAGE, "NORMALIZE must be Boolean."
        ENDIF

        NORMALIZE = LONG(NORMALIZE)

    ENDIF ELSE BEGIN
        NORMALIZE = 0L
    ENDELSE

    ; Construct the filename with details about the data set.
    OUTWTYPE = STRCOMPRESS(STRUCT.WTYPE, /REMOVE_ALL)
    DATE = DATE_STRING()
    OUTRES = LONG(1000L * STRUCT.RES)
    RESDSC = 'RES_' + STRCOMPRESS(STRING(OUTRES), /REMOVE_ALL) + 'm'
    OUTRANGE = MINMAX(STRUCT.RHO_KM_VALS)
    LNRNGA = STRCOMPRESS(STRING(LONG(OUTRANGE[0])), /REMOVE_ALL)
    LNRNGB = STRCOMPRESS(STRING(LONG(OUTRANGE[1])), /REMOVE_ALL)
    RNGDSC = 'RNG_' + LNRNGA + '_' + LNRNGB + 'km'
    OUTSTR = DIR + DATE + '_' + RESDSC + '_' + RNGDSC + '_' + OUTWTYPE

    IF TAG_EXIST(STRUCT, 'P_NORM_FWD_VALS') THEN OUTSTR += '_FWD'
    IF TAG_EXIST(STRUCT, 'POWER_VALS') THEN OUTSTR += '_INT'
    IF TAG_EXIST(STRUCT, 'POWER_FFT_VALS') THEN OUTSTR += '_FFT'

    CASE LAGRANGE OF
        0L: OUTSTR += '_NORMAL'
        1L: OUTSTR += '_LAGRANGE'
        2L: OUTSTR += '_FAST'
    ENDCASE

    IF NORMALIZE THEN OUTSTR +='_NORM'
    OUTSTR += '_' + REV + '.sav'
    RETURN, OUTSTR
END
