cmake_minimum_required(VERSION 3.20)
include_guard()
include(FetchContent)

if(MSVC)
    set(CERF_CPP ON)
else()
    set(CERF_CPP OFF)
endif()
set(LIB_MAN OFF)

fetchcontent_declare(
    cerf_checkout
    GIT_REPOSITORY https://jugit.fz-juelich.de/mlz/libcerf.git
    GIT_TAG
        4f2c6384e12d51ea86302bab0c0cb53b2e948c19 # the commit corresponding to v1.17
)

fetchcontent_makeavailable(cerf_checkout)

if(cerf_checkout_POPULATED)
    set(cerf_INCLUDE_DIRECTORY ${cerf_checkout_SOURCE_DIR}/lib)
else()
    message(FATAL_ERROR "cerf is not available")
endif()
