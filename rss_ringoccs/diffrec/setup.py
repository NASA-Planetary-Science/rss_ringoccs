def configuration(parent_package=None, top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration(package_name=None,
                           package_path=None,
                           parent_name=parent_package,
                           top_path=top_path)
    config.add_extension(
        '_special_functions', [
            "src/__bessel_J0.c",
            "src/__bessel_I0.c",
            "src/__fraunhofer_diffraction_double_slit.c",
            "src/__fraunhofer_diffraction_single_slit.c",
            "src/__fresnel_diffraction_inverted_square_well.c",
            "src/__fresnel_diffraction_left_straightedge.c",
            "src/__fresnel_diffraction_right_straightedge.c",
            "src/__fresnel_diffraction_square_well_phase.c",
            "src/__fresnel_diffraction_square_well.c",
            "src/__fresnel_integral_cosine.c",
            "src/__fresnel_integral_cosine_heald.c",
            "src/__fresnel_integral_cosine_while.c",
            "src/__fresnel_integral_sine.c",
            "src/__fresnel_integral_sine_heald.c",
            "src/__fresnel_integral_sine_while.c",
            "src/__fresnel_integral.c",
            "src/__fresnel_kernel.c",
            "src/__get_array.c",
            "src/__math_function_lambertw.c",
            "src/__math_function_legendre.c",
            "src/__math_function_max.c",
            "src/__math_function_min.c",
            "src/__math_function_norm_eq_width.c",
            "src/__math_function_resolution_inverse.c",
            "src/__math_function_sinc.c",
            "src/_special_functions.c",
            "src/__where_greater.c",
            "src/__where_lesser.c"
        ]
    )
    config.add_extension('_window_functions', [
            "src/_window_functions.c",
            "src/__bessel_I0.c",
            "src/__window_function_squared_cosine.c",
            "src/__window_function_kaiser_bessel.c",
            "src/__window_function_normalization.c",
            "src/__window_function_modified_kaiser_bessel.c",
            "src/__window_function_rectangular.c"
        ]
    )
    config.add_extension('_diffraction_functions', [
            "src/_diffraction_functions.c",
            "src/__diffraction_functions.c",
            "src/__fresnel_kernel.c",
            "src/__bessel_I0.c",
            "src/__math_function_legendre.c",
            "src/__fresnel_transform.c",
            "src/__fresnel_transform_ellipse.c",
            "src/__fresnel_transform_legendre.c",
            "src/__fresnel_transform_newton.c",
            "src/__window_function_squared_cosine.c",
            "src/__window_function_kaiser_bessel.c",
            "src/__window_function_modified_kaiser_bessel.c",
            "src/__window_function_rectangular.c"
        ]
    )

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)
