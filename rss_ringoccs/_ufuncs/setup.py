def configuration(parent_package=None, top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration(package_name='rss_ringoccs',
                           package_path=None,
                           parent_name=parent_package,
                           top_path=top_path)
    config.add_extension('_special_functions', ['_special_functions.c'])
    config.add_extension('_fresnel_integrals', ['_fresnel_integrals.c'])

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)

