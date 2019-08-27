def configuration(parent_package=None, top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration(package_name=None,
                           package_path=None,
                           parent_name=parent_package,
                           top_path=top_path)
    config.add_extension('_special_functions',
                         ['src/_special_functions.c'])
    config.add_extension('_window_functions',
                         ['src/_window_functions.c'])
    config.add_extension('_diffraction_functions',
                         ['src/_diffraction_functions.c'])

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)
