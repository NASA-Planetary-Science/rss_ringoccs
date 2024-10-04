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
;       Checks to make sure the geometry variables are correct.                ;
;   Arguments:                                                                 ;
;       INSTRUCT (structure):                                                  ;
;           Structure with the reconstruction data.                            ;
;   Output:                                                                    ;
;       OCCSTRUCT (string):                                                    ;
;           The corrected INSTRUCT.                                            ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018/04/04                                                         ;
;------------------------------------------------------------------------------;
FUNCTION GET_OCC_TYPE, INSTRUCT

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    IF IDL_TYPE(INSTRUCT) NE IDLTYPE({A:'Struc'}) THEN BEGIN
        MESSAGE, "Input must be a structure."
    ENDIF

    DRHO = MINMAX(INSTRUCT.RHO_DOT_KMS_VALS)
    DX = INSTRUCT.RHO_KM_VALS[1] - INSTRUCT.RHO_KM_VALS[0]

    ;Check to make sure there is no ambiguity about the type of occultation.
    IF (DRHO[0] LT 0) AND (DRHO[1] GT 0) THEN BEGIN
        MESSAGE, "ERROR: dr/dt is positive and negative. This is a chord occ."
    ENDIF

    IF (DX GT 0) AND (DRHO[1] LT 0) THEN BEGIN
        INSTRUCT.RHO_DOT_KMS_VALS = ABS(INSTRUCT.RHO_DOT_KMS_VALS)
    ENDIF

    IF DX LT 0 THEN BEGIN
        INSTRUCT.RHO_KM_VALS = REVERSE(INSTRUCT.RHO_KM_VALS)
        INSTRUCT.PHASE_RAD_VALS = REVERSE(INSTRUCT.PHASE_RAD_VALS)
        INSTRUCT.RAW_POWER_VALS = REVERSE(INSTRUCT.RAW_POWER_VALS)
        INSTRUCT.PHI_RAD_VALS = REVERSE(INSTRUCT.PHI_RAD_VALS)
        INSTRUCT.B_RAD_VALS = REVERSE(INSTRUCT.B_RAD_VALS)
        INSTRUCT.F_SKY_HZ_VALS = REVERSE(INSTRUCT.F_SKY_HZ_VALS)
        INSTRUCT.D_KM_VALS = REVERSE(INSTRUCT.D_KM_VALS)
        INSTRUCT.RHO_DOT_KMS_VALS = ABS(REVERSE(INSTRUCT.RHO_DOT_KMS_VALS))
    ENDIF

    RETURN, INSTRUCT
END
