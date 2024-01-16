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
;       Computes the product of the window function with the Fresnel kernel.   ;
;   Arguments:                                                                 ;
;       WINDOW_ARRAY (real, array-like):                                       ;
;           The window function.                                               ;
;       PSI_ARRAY (real, array-like):                                          ;
;           The Fresnel kernel.                                                ;
;   Output:                                                                    ;
;       WKERNEL (complex, array-like):                                         ;
;           The windowed Fresnel kernel.                                       ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2018/04/12                                                         ;
;------------------------------------------------------------------------------;

; Computes the product of the window function and the Fresnel kernel.
FUNCTION WINDOW_KERNEL, WINDOW_ARRAY, PSI_ARRAY

    ; Error checking codes.
    ON_ERROR, 2

    ; The complex Fresnel kernel is e^{i psi}. Compute i psi.
    IPSI_ARRAY = COMPLEX(0.0, PSI_ARRAY)

    ; Compute the point-wise product of the window function and the kernel.
    WKERNEL = WINDOW_ARRAY * EXP(IPSI_ARRAY)
    RETURN, WKERNEL
END
