def configuration(parent_package=None, top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration(package_name=None,
                           package_path=None,
                           parent_name=parent_package,
                           top_path=top_path)
    config.add_extension(
        '_special_functions', [
            "src/__sinc.c",
            "src/__fraunhofer_diffraction_double_slit.c",
            "src/__fraunhofer_diffraction_single_slit.c",
            "src/__fresnel_diffraction_inverted_square_well.c",
            "src/__fresnel_diffraction_left_straightedge.c",
            "src/__fresnel_diffraction_right_straightedge.c",
            "src/__fresnel_diffraction_square_well_phase.c",
            "src/__fresnel_diffraction_square_well.c",
            "src/__fresnel_integral_cosine.c",
            "src/__fresnel_integral_sine.c",
            "src/__fresnel_integral.c",
            "src/__fresnel_kernel.c",
            "src/__math_function_bessel.c",
            "src/__math_function_lambertw.c",
            "src/__math_function_legendre.c",
            "src/__math_function_max.c",
            "src/__math_function_min.c",
            "src/__math_function_norm_eq_width.c",
            "src/__math_function_resolution_inverse.c",
            "src/_special_functions.c"
        ]
    )
    config.add_extension('_window_functions', [
            "src/_window_functions.c",
            "src/__window_function_squared_cosine.c",
            "src/__window_function_kaiser_bessel.c",
            "src/__window_function_modified_kaiser_bessel.c",
            "src/__window_function_rectangular.c"
        ]
    )
    config.add_extension('_diffraction_functions', [
            "src/_diffraction_functions.c",
            "src/__math_function_legendre.c",
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
