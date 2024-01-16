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
;       Computes the Lambert W function.                                       ;
;   Arguments:                                                                 ;
;       XIN (real, array-like):                                                ;
;           The independent variable.                                          ;
;   Output:                                                                    ;
;       LWX (real, array-like):                                                ;
;           The Lambert W function of XIN.                                     ;
;   References:                                                                ;
;       1.) Gonnet, Hare, Jeffrey, Knuth, 1996                                 ;
;           On the Lambert W function.                                         ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018                                                               ;
;------------------------------------------------------------------------------;

; Function for computing the Lambert W function.
FUNCTION LAMBERTW, XIN

    ; Error checking code.
    ON_ERROR, 2

    ; The inputs should be real valued.
    IF (CHECK_REAL(XIN) EQ 0) THEN BEGIN
        MESSAGE, "Independent variable should be real-valued."
    ENDIF

    ; A good initial guess is x / exp(1). Compute this.
    X0 = XIN * 0.36787944

    ; Desired tolerance, single precision.
    TOLER = 1.0E-7

    ; The exponential constant, exp(1).
    EULER_E = 2.7182818

    ; For x > e, a better guess is log(x) - log(log(x)). Check if we need this.
    BIG_IND = WHERE(XIN GT EULER_E)

    ; For any points with real part larger than e use x0 = log(x) - log(log(x)).
    IF (MIN(BIG_IND) NE -1) THEN BEGIN
        ; We can avoid a call to log by precomputing log(x).
        LOG_XIN_BIG = ALOG(X0[BIG_IND])

        ; log(x) - log(log(x)) now needs only one more log call.
        X0[BIG_IND] = LOG_XIN_BIG - ALOG(LOG_XIN_BIG)
    ENDIF

    ; Use Halley's method to finish the computation.
    RETURN, LAMBERTW_HALLEY(XIN, X0, TOLER)
END
