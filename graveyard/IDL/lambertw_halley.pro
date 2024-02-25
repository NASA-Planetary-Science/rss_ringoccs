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
;       Computes the Lambert W function via Halley's method.                   ;
;   Arguments:                                                                 ;
;       XIN (real, array-like):                                                ;
;           The independent variable.                                          ;
;       X0 (real, array-like):                                                 ;
;           The initial guess for Halley's method.                             ;
;       TOL (real):                                                            ;
;           The allowed tolerance in Halley's method.                          ;
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

; Computes the Lambert W function using Halley's method.
FUNCTION LAMBERTW_HALLEY, XIN, X0, TOL

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; The inputs should be real valued.
    IF (CHECK_REAL(XIN) EQ 0) THEN BEGIN
        MESSAGE, "Independent variable should be real-valued."
    ENDIF

    ; Same check for the initial guess term.
    IF (CHECK_REAL(X0) EQ 0) THEN BEGIN
        MESSAGE, "Initial guess should be real-valued."
    ENDIF

    ; The tolerance should be a strictly-positive real number.
    IF (CHECK_POSITIVE_REAL_NUMBER(TOL) EQ 0) THEN BEGIN
        MESSAGE, "Tolerance should be a positive real number."
    ENDIF

    ; The Halley term involves the derivatives of x*exp(x). Compute these terms.
    EXP_X0 = EXP(X0)
    S = X0 + 1.0
    T = X0*EXP_X0 - XIN

    ; The Incremental factor for Halley's method. X_{n+1} = X_{n} - DX.
    DX = T / (EXP_X0*S - 0.5*(S + 1.0)*T / S)

    ; Iteratively perform Halley's method. 10 iterations is more than enough.
    FOR N = 0, 10 DO BEGIN

        ; If all entries of the array have converged to the desire tolerance,
        ; we can break out of the main for loop.
        IF (MAX(ABS(DX)) LT TOL) THEN BEGIN
            BREAK
        ENDIF

        ; Otherwise, update the Halley iterate to the next term.
        X0 = X0 - DX

        ; Compute the derivative factors at the new guess.
        EXP_X0 = EXP(X0)
        S = X0 + 1.0
        T = X0*EXP_X0 - XIN

        ; And lastly, update the incremental factor.
        DX = T / (EXP_X0*S - 0.5*(S + 1.0)*T / S)
    ENDFOR

    ; The Halley increment is X_{n+1} = X_{n} - DX. Return this.
    RETURN, X0 - DX
END
