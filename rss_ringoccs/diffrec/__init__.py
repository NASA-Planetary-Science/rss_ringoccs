"""
Copyright (C) 2018 Team Cassini

This program is free software: you  can redistribute it and/or modify it
under the  terms of the GNU  General Public License as  published by the
Free Software Foundation,  either version 3 of the License,  or (at your
option) any later version.

This  program  is distributed  in  the  hope  that  it will  be  useful,
but WITHOUT  ANY  WARRANTY;  without   even  the  implied  warranty  of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy  of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.

This program is part of the rss_ringoccs repository hosted at
https://github.com/NASA-Planetary-Science/rss_ringoccs and developed
with the financial support of NASA's Cassini Mission to Saturn.

    Subpackage Name:
        diffrec
    Purpose:
        Provide functions and classes that aid in the process of
        Diffraction Correction / Fresnel Inversion. Additional
        functions for the purpose of forward modelling of
        reconstructed data and diffraction modelling are included.
        Special mathematical functions and solutions to common
        problems in diffraction theory are also included.
    Sub-Modules:
        advanced_tools:
            This submodule is good for modeling the geometry of a
            given occultation, and for comparing your results to
            the results obtained by others (Ex. Results on the PDS).
            This submodule contains the following:
                compare_tau:
                    Class
                    Used for running diffraction correction on a
                    given set of diffracted data and then comparing
                    the outcome to a given set of reconstructed data.
                find_optimal_resolution:
                    Class
                    Given a set of data and a reconstruction, this
                    class will run diffraction correction over a set
                    of resolutions. The output contains the L_2 and
                    L_infinity difference of the two reconstructions
                    as a function of resolution.
                delta_impulse_diffraction:
                    Class
                    Given a set of geometry data, this class create
                    modeled data of the solution of diffraction
                    through a Dirac-Delta impulse function.
                    Reconstruction is then perform on the mock-data.
                    This tool is good for modeling the problem and
                    determining resolution constraints based on the
                    geometry available.
        diffraction_correction:
            This is the main sub-module in the entire subpackage.
            Given a set of diffracted data and a requested resolution
            (in kilometers), diffraction corrections will be
            performed to produce a diffraction corrected profile.
            This submodule comtains the following:
                DiffractionCorrection:
                    Class
                    Given a requested resolution and an instance of
                    the NormDiff class (See Calibration subpackage),
                    this produces a diffraction corrected profile.
        Special Functions:
            fresnel_sin.........The Fresnel sine integral.
            fresnel_cos.........The Fresnel cosine integral.
            sq_well_solve.......Diffraction pattern through square well.
            compute_norm_eq.....Computes the normalized equivalent width.
            resolution_inverse..Computes the inverse of the function
                                y = x/(exp(-x)+x-1)
            fresnel_scale.......Compute the Fresnel scale.
        window_functions:
            rect................Rectangular window.
            coss................Squared cossine window.
            kb20................Kaiser-Bessel 2.0 window.
            kb25................Kaiser-Bessel 2.5 window.
            kb35................Kaiser-Bessel 3.5 window.
            kbmd20..............Modified Kaiser-Bessel 2.0 window.
            kbmd25..............Modified Kaiser-Bessel 2.5 window.
            get_range_actual....Given an array of numbers (usually the
                                radial range of the data), a range
                                request, and a window width, compute the
                                allowed range of processing.
"""

# Import the main classes used in diffraction reconstruction.
from .diffraction_correction import DiffractionCorrection
