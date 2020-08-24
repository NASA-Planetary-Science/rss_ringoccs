def configuration(parent_package=None, top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration(package_name=None,
                           package_path=None,
                           parent_name=parent_package,
                           top_path=top_path)
    config.add_extension(
        '_special_functions', [
            "src/_bessel_I0.c",
            "src/_bessel_J0.c",
            "src/_diffraction_functions.c",
            "src/_fraunhofer_diffraction.c",
            "src/_frequency_to_wavelength.c",
            "src/_fresnel_diffraction_gap.c",
            "src/_fresnel_diffraction_left_straightedge.c",
            "src/_fresnel_diffraction_right_straightedge.c",
            "src/_fresnel_diffraction_ringlet_phase.c",
            "src/_fresnel_diffraction_ringlet.c",
            "src/_fresnel_diffraction_square_wave.c",
            "src/_fresnel_integral_cosine.c",
            "src/_fresnel_integral_sine.c",
            "src/_fresnel_integral.c",
            "src/_fresnel_kernel.c",
            "src/_fresnel_scale.c",
            "src/_fresnel_transform.c",
            "src/_fresnel_transform_ellipse.c",
            "src/_fresnel_transform_interpolate.c",
            "src/_fresnel_transform_legendre_even.c",
            "src/_fresnel_transform_legendre_odd.c",
            "src/_fresnel_transform_newton.c",
            "src/_fresnel_transform_perturbednewton.c",
            "src/_lambertw.c",
            "src/_legendre.c",
            "src/_max.c",
            "src/_min.c",
            "src/_norm_eq_width.c",
            "src/_resolution_inverse.c",
            "src/_sinc.c",
            "src/_wavelength_to_wavenumber.c",
            "src/_where.c",
            "src/_window_function_squared_cosine.c",
            "src/_window_function_kaiser_bessel.c",
            "src/_window_function_modified_kaiser_bessel.c",
            "src/_window_function_rectangular.c",
            "src/_window_function_normalization.c",
            "src/special_functions.c",
        ]
    )

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)
