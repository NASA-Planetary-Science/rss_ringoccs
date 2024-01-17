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
;       Interpolates the Fresnel kernel against a data set.                    ;
;   Arguments:                                                                 ;
;       XVALS (real, array-like):                                              ;
;           The data points for psi.                                           ;
;       XINTERP (real, array-like):                                            ;
;           The desired points to interpolate psi to.                          ;
;       FRESNEL_PSI (real, array-like):                                        ;
;           The unexponentiated Fresnel kernel.                                ;
;   Outputs:                                                                   ;
;       FRESNEL_KERNEL (complex, array-like):                                  ;
;           The complex interpolated Fresnel kernel.                           ;
;------------------------------------------------------------------------------;
;   Author: Ryan Maguire                                                       ;
;   Date:   2017                                                               ;
;------------------------------------------------------------------------------;

; Interpolates the Fresnel kernel and returns the complex form.
FUNCTION FRESNEL_KERNEL_INTERPOLATE, XVALS, XINTERP, FRESNEL_PSI

    ; Interpolate the data against the desired set of points.
    PSI_INTERP = INTERPOL(FRESNEL_PSI, XINTERP, XVALS)

    ; Compute the complex Fresnel kernel and return.
    FRESNEL_KERNEL = EXP(COMPLEX(0.0, FRESNEL_PSI))
    RETURN, FRESNEL_KERNEL
END
