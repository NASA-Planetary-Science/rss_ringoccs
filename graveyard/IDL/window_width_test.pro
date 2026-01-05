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
MAX_ABS_DIFFERENCE_KM = FLTARR(20)
DISPLACEMENT = MAX_ABS_DIFFERENCE_KM / 1.d0

FOR i = 0, N_ELEMENTS(MAX_ABS_DIFFERENCE_KM) - 1 DO BEGIN
  WINDOW_KM_REQUEST = 30.d0 + DISPLACEMENT[i]
  DIFFRACTION_CORRECTION_WITH_WINDOW_DIFFERENCE, WINDOW_KM_REQUEST, DIFFERENCE
  MAX_ABS_DIFFERENCE_KM[i] = DIFFERENCE
ENDFOR

END
