cmake_minimum_required(VERSION 3.20)
include_guard()
include(FetchContent)
set(FETCHCONTENT_QUIET OFF)

set(LIBTMPL_BUILD_EXAMPLES OFF)
set(LIBTMPL_BUILD_TESTS OFF)
set(LIBTMPL_USE_ASM ON)
set(LIBTMPL_USE_MATH ON)
set(LIBTMPL_STDVER "-std=c11")
if(WIN32)
    set(LIBTMPL_USE_FASM OFF)
    set(LIBTMPL_USE_OMP OFF)
elseif(APPLE)
    set(LIBTMPL_USE_FASM OFF)
    set(LIBTMPL_USE_OMP ON)
else()
    set(LIBTMPL_USE_FASM ON)
    set(LIBTMPL_USE_OMP ON)
endif()

fetchcontent_declare(
    tmpl_checkout
    GIT_REPOSITORY https://github.com/tcumby/libtmpl.git
    GIT_TAG add-cmake-build-system-redo
    SOURCE_DIR libtmpl
)

fetchcontent_makeavailable(tmpl_checkout)
