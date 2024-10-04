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
;       This function computes the possible allowed range for processing.      ;
;   Arguments:                                                                 ;
;       RHO (array):                                                           ;
;           The radial range of the data.                                      ;
;       RANGE (array):                                                         ;
;           The requested start and end points for processing, [a, b].         ;
;       W_KM_VALS (array):                                                     ;
;           The window size as a function of rho.                              ;
;   Keywords:                                                                  ;
;       N_USED (int):                                                          ;
;           The number of points in the processed range.                       ;
;       FINISH (int):                                                          ;
;           The final index of the processing range.                           ;
;   Output:                                                                    ;
;       START (int):                                                           ;
;           The starting index for the range.                                  ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018/04/10                                                         ;
;------------------------------------------------------------------------------;
FUNCTION GET_RANGE_ACTUAL, RHO, RANGE, W_KM_VALS, $
                           N_USED = N_USED, FINISH = FINISH

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    IF ~CHECK_REAL(RHO) THEN BEGIN
        MESSAGE, 'RHO must be array of real numbers.'
    ENDIF

    IF MIN(RHO) LT 0.d0 THEN BEGIN
        MESSAGE, 'RHO must be positive.'
    ENDIF

    IF N_ELEMENTS(RANGE) NE 2L THEN BEGIN
        MESSAGE, 'RANGE must have format RANGE = [a, b].'
    ENDIF

    IF ~CHECK_POS_REAL(MIN(RANGE)) THEN BEGIN
        MESSAGE, 'RANGE must be positive.'
    ENDIF

    IF ~CHECK_POS_REAL(MAX(RANGE)) THEN BEGIN
        MESSAGE, 'RANGE must be positive.'
    ENDIF

    IF ~CHECK_REAL(W_KM_VALS) THEN BEGIN
        MESSAGE, 'W_KM_VALS must be real.'
    ENDIF

    IF MIN(W_KM_VALS) LT 0.0 THEN BEGIN
        MESSAGE, 'W_KM_VALS must be positive.'
    ENDIF

    IF MIN(RANGE) GT MAX(RHO) THEN BEGIN
        MESSAGE, 'Data not available: MIN(request) > MAX(actual).'
    ENDIF

    IF MAX(RANGE) LT MIN(RHO) THEN BEGIN
        MESSAGE, 'Data not available: MAX(request) < MIN(actual).'
    ENDIF

    W_MAX = MAX(W_KM_VALS)

    RHO_MIN_LIM = MIN(RHO) + CEIL(0.5 * W_MAX)
    RHO_MAX_LIM = MAX(RHO) - CEIL(0.5 * W_MAX)

    RHO_START = RHO[MIN(WHERE(RHO GE MIN(RANGE)))]
    RHO_END = RHO[MAX(WHERE(RHO LE MAX(RANGE)))]

    RHO_MIN = MAX([RHO_MIN_LIM, RHO_START])
    RHO_MAX = MIN([RHO_MAX_LIM, RHO_END])

    START = LONG(MIN(WHERE(RHO GE RHO_MIN)))
    FINISH = LONG(MAX(WHERE(RHO LE RHO_MAX)))
    N_USED = FINISH - START

    RETURN, START
    IF KEYWORD_SET(N_USED) THEN RETURN, N_USED
    IF KEYWORD_SET(FINISH) THEN RETURN, FINISH
END
