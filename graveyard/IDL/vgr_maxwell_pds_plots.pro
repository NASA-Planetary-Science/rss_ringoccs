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
;       Compares our inversion of the Maxwell ringlet from Voyager to the PDS. ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Tells the compiler that integers should be 32 bits, not 16.
COMPILE_OPT IDL2

; The data file.
RAWFILE = 'VG_2803_peer_review/S_RINGS/RAWDATA/RS1R1BDP.DAT'

; Number of variables in the file.
NUMBER_OF_COLUMNS = 16

; Number of rows.
NUMBER_OF_POINTS = 375050L

; Distance between points in the radial direction.
SAMPLE_SPACING = 0.20

; Radii of the points, in kilometers.
RADIUS_KM = 70000.0 + FINDGEN(NUMBER_OF_COLUMNS) * SAMPLE_SPACING

; Empty array for the data to be stored in.
DATA_RAW = FLTARR(NUMBER_OF_COLUMNS, NUMBER_OF_POINTS)

; Read in the data file and load it into the data buffer.
CLOSE, 1
OPENR, 1, RAWFILE
READU, 1, RAWFILE
CLOSE, 1

; The variables we care about.
TAU_RAW = DATA_RAW[0, *]
TAU_PDS_RAW = DATA_RAW[4, *]

; Ring opening angle, in radians.
B_ANGLE = 84.07460 * !PI / 180.0

; Normalization factor for the power.
MU = COS(B_ANGLE)

; The power for the PDS and our profile.
POWER_RAW = EXP(-TAU_RAW / MU)
POWER_PDS = EXP(-TAU_PDS_RAW / MU)

; Make some plots.
NSMOOTH = 5
PLOT, RADIUS_KM, SMOOTH(POWER_RAW, NSMOOTH)
OPLOT, RADIUS_KM, SMOOTH(POWER_PDS, NSMOOTH), COLOR = CGCOLOR('red')

END
