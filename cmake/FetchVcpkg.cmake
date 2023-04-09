cmake_minimum_required(VERSION 3.20)
include_guard()
include(FetchContent)
set(FETCHCONTENT_QUIET OFF)
fetchcontent_declare(
    vcpkg_checkout
    GIT_REPOSITORY https://github.com/microsoft/vcpkg.git
    GIT_TAG
        c64c0fdac572ca43ea5ae018fc408ddced50d5b1 # 2022.02.02
)

fetchcontent_makeavailable(vcpkg_checkout)

if(vcpkg_checkout_POPULATED)
    set(vcpkg_toolchain_file
        "${vcpkg_checkout_SOURCE_DIR}/scripts/buildsystems/vcpkg.cmake"
    )
else()
    message(FATAL_ERROR "Failed to checkout vcpkg.")
endif()
