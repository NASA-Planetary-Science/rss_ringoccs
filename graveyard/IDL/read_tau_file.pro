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
;       Reads in the data from a Tau file.                                     ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Extracts the data from a Tau file.
PRO READ_TAU_FILE, TAU_FILE, TAU_STRUCT

    ; Tells the compiler that integers should be 32 bits, not 16.
    COMPILE_OPT IDL2

    ; Error checking code.
    ON_ERROR, 2

    ; Read in the data file into an array.
    TAU_IN = RDCOL(TAU_FILE, 1, 99999999, [1:11])

    ; Label the columns for the array in the struct to return.
    TAU_STRUCT = {                  $
        RHO_KM:TAU_IN[0, *],        $
        RHO_COR_KM:TAU_IN[1, *],    $
        PHI_RL_DEG:TAU_IN[2, *],    $
        PHI_ORA_DEG:TAU_IN[3, *],   $
        TAU_RECON:TAU_IN[4, *],     $
        PHASE_DEG:TAU_IN[5, *],     $
        T_OET_SPM:TAU_IN[7, *],     $
        T_RET_SPM:TAU_IN[8, *],     $
        T_SET_SPM:TAU_IN[9, *],     $
        B_DEG:TAU_IN[10, *]         $
    }

END
