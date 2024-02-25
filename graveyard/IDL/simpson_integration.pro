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
;       Numerical integration using Simpson's method.                          ;
;   Arguments:                                                                 ;
;       XIN (real, array-like):                                                ;
;           The (evenly spaced) independent variable.                          ;
;       YIN (real, array-like):                                                ;
;           The data corresponding to XIN.                                     ;
;   Output:                                                                    ;
;       SUM (real):                                                            ;
;           The numerical integral of YIN against XIN.                         ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Function for computing the numerical integral of a data set.
FUNCTION SIMPSON_INTEGRATION, XIN, YIN

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; The data sets should be real-valued.
    IF (CHECK_REAL(XIN) EQ 0) OR (CHECK_REAL(YIN) EQ 0) THEN BEGIN
        MESSAGE, "Input data should be real-valued."
    ENDIF

    ; Number of elements in the data set.
    X_LENGTH = N_ELEMENTS(XIN)

    ; Center of the data. We use this to split the data set in to even and
    ; odd points to perform the summation.
    MID_POINT = (X_LENGTH - 2) / 2

    ; This assumes evenly spaced samples.
    DX = (XIN[-1] - XIN[0]) / (FLOAT(X_LENGTH) - 1.0)

    ; Multiplicative factor for the Simpson sum.
    FACTOR = DX / 3.0

    ; Two cases: Even number of points and odd. The weights are handled
    ; differently in each case. First, the even case.
    IF IS_EVEN(X_LENGTH) THEN BEGIN

        ; Make the data into a 2x(N/2) array.
        ZIN = REFORM(YIN[1:X_LENGTH - 2], 2, MID_POINT)

        ; Sum down the columns. The zeroth entry is the sum of the even points
        ; and the first entry is the sum of the odd points.
        ZSUM = TOTAL(ZIN, 2)

        ; The Simpson sum can now be performed by multipying by the weights.
        SUM = YIN[0] + 2.0*ZSUM[0] + 4.0*ZSUM[1] + YIN[X_LENGTH - 1]

    ; Second case, odd number of points.
    ENDIF ELSE BEGIN

        ; Same split as before, but neglect the zeroth data point. We are then
        ; looking at an even number of points and use the previous sum.
        ZIN = REFORM(YIN[2:X_LENGTH - 2], 2, MID_POINT)

        ; Sum down the columns again.
        ZSUM = TOTAL(ZIN, 2)

        ; Use the triangle approximation for the first two data points, and
        ; then Simpson's formula for the rest.
        SUM = YIN[0] + 4.0*YIN[1] + 2.0*ZSUM[0] + 4.0*ZSUM[1] + YIN[X_LENGTH-1]
    ENDELSE

    ; The numerical integral is the weighted sum. Multiply and return.
    RETURN, SUM * FACTOR
END
