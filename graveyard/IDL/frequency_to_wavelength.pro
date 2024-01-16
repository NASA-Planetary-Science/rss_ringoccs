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
;       Converts frequency in Hertz to wavelength in kilometers.               ;
;   Arguments:                                                                 ;
;       FREQUENCY (real, array-like):                                          ;
;           An array of positive inputs, the frequency of the light.           ;
;   Output:                                                                    ;
;       WAVELENGTH (real, array-like):                                         ;
;           The corresponding wavelength, in kilometers.                       ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018/04/18                                                         ;
;------------------------------------------------------------------------------;

; Function for converting frequency to wavelength for light in a vacuum.
FUNCTION FREQUENCY_TO_WAVELENGTH, FREQUENCY

    ; Error checking code.
    ON_ERROR, 2

    ; Check for invalid inputs. Should have at least one element.
    IF N_ELEMENTS(FREQUENCY) EQ 0 THEN BEGIN
        MESSAGE, 'Input is undefined.'
    ENDIF

    ; Frequency should be real-valued (and in Hertz).
    IF CHECK_REAL(FREQUENCY) EQ 0 THEN BEGIN
        MESSAGE, 'Input must be real.'
    ENDIF

    ; Negative and zero-valued frequencies are undefined.
    IF MIN(FREQUENCY) LE 0.0 THEN BEGIN
        MESSAGE, 'Input must be positive.'
    ENDIF

    ; The speed of light in kilometers per second.
    SPEED_OF_LIGHT = 299792.458

    ; The wavelength is inversely related to the frequency. Compute this.
    RETURN, SPEED_OF_LIGHT / FREQUENCY
END
