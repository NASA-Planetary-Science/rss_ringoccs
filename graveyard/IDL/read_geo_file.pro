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
;       Reads in the data from a Geo file.                                     ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Extracts the data from a Geo file.
PRO READ_GEO_FILE, GEO_FILE, GEO_STRUCT

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; Read in the data file into an array.
    GEO_IN = rdcol(GEO_FILE, 1, 99999999, [1:18])

    ; Label the columns for the array in the struct to return.
    GEO_STRUCT = {                      $
        T_OET_SPM:GEO_IN[0, *],         $
        T_RET_SPM:GEO_IN[1, *],         $
        T_SET_SPM:GEO_IN[2, *],         $
        RHO_KM:GEO_IN[3, *],            $
        PHI_RL_DEG:GEO_IN[4, *],        $
        PHI_ORA_DEG:GEO_IN[5, *],       $
        B_DEG:GEO_IN[6, *],             $
        D_KM:GEO_IN[7, *],              $
        RHO_DOT_KMS:GEO_IN[8, *],       $
        PHI_RL_DOT_KMS:GEO_IN[9, *],    $
        F_KM:GEO_IN[10, *],             $
        R_IMP_KM:GEO_IN[11, *],         $
        R_SC_X_KM:GEO_IN[12, *],        $
        R_SC_Y_KM:GEO_IN[13, *],        $
        R_SC_Z_KM:GEO_IN[14, *],        $
        R_SC_DOT_X_KMS:GEO_IN[15, *],   $
        R_SC_DOT_Y_KMS:GEO_IN[16, *],   $
        R_SC_DOT_Z_KMS:GEO_IN[17, *]    $
    }

END
