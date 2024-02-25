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
;       Computes the Fresnel scale from real-valued inputs.                    ;
;   Arguments:                                                                 ;
;       LAMBDA (real, array-like):                                             ;
;           The wavelength of the incident ray.                                ;
;       D (real, array-like):                                                  ;
;           The distance from the observer to the ring point.                  ;
;       PHI (real, array-like):                                                ;
;           The ring azimuth angle.                                            ;
;       B (real, array-like):                                                  ;
;           The ring opening angle.                                            ;
;   Keywords:                                                                  ;
;       DEG (Boolean):                                                         ;
;           Boolean indicating if angles are in degrees or radians.            ;
;   Output:                                                                    ;
;       FRES (real, array-like):                                               ;
;           The Fresnel scale.                                                 ;
;   References:                                                                ;
;       1.) Marouf, Tyler, Rosen, 1986                                         ;
;           Profiling Saturn's Rings by Radio Occultation.                     ;
;       2.) Goodman, 1969                                                      ;
;           Introduction to Fourier Optics.                                    ;
;       3.) Carroll, Ostie, 2007                                               ;
;           An Introduction to Modern Astrophysics.                            ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018/04/04                                                         ;
;------------------------------------------------------------------------------;

; Computes the Fresnel scale from real-valued inputs.
FUNCTION FRESNEL_SCALE, LAMBDA, D, PHI, B, DEG = DEG

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR,2

    ; The wavelength should be real-valued (and positive).
    IF (CHECK_REAL(LAMBDA) EQ 0L) THEN BEGIN
        MESSAGE, 'LAMBDA must be real'
    ENDIF

    ; Distance from the spacecraft to the rings. Must be real.
    IF (CHECK_REAL(D) EQ 0L) THEN BEGIN
        MESSAGE, 'D must be real'
    ENDIF

    ; The azimuthal angle of the point in the plane. Must be real.
    IF (CHECK_REAL(PHI) EQ 0L) THEN BEGIN
        MESSAGE, 'PHI must be real'
    ENDIF

    ; The ring opening angle with respect to Earth. Must be real.
    IF (CHECK_REAL(B) EQ 0L) THEN BEGIN
        MESSAGE, 'B must be real'
    ENDIF

    ; The input is allowed to be degrees or radians.
    IF KEYWORD_SET(DEG) THEN BEGIN

        ; Ensure the user didn't give an invalid keyword for the DEG variable.
        IF (CHECK_BOOLE(DEG) EQ 0) THEN BEGIN
            MESSAGE, "DEG must be Boolean"
        ENDIF

        ; If degrees were requested, use the corresponding trig functions.
        IF DEG THEN BEGIN
            COS_B = COSD(B)
            SIN_B = SIND(B)
            SIN_P = SIND(PHI)

        ; Otherwise perform the same calculation in radians.
        ENDIF ELSE BEGIN
            COS_B = COS(B)
            SIN_B = SIN(B)
            SIN_P = SIN(PHI)
        ENDELSE

    ; Default is radians.
    ENDIF ELSE BEGIN
        COS_B = COS(B)
        SIN_B = SIN(B)
        SIN_P = SIN(PHI)
    ENDELSE

    ; The Fresnel scale makes frequent use of the square of the previous values.
    COS_B_SQ = COS_B * COS_B
    SIN_B_SQ = SIN_B * SIN_B
    SIN_P_SQ = SIN_P * SIN_P

    ; Factors for the square of the unnormalized Fresnel scale.
    NUM_FACTOR = 1.0 - COS_B_SQ * SIN_P_SQ
    DEN_FACTOR = 2.0 * SIN_B_SQ

    ; The Fresnel scale normalizes by the wavelength and distance from the
    ; observer to the point of interest. Compute this and take square roots.
    FRES = SQRT(LAMBDA * D * NUM_FACTOR / DEN_FACTOR)
    RETURN, FRES
END
