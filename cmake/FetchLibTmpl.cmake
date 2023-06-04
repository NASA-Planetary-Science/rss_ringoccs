cmake_minimum_required(VERSION 3.20)
include_guard()
include(FetchContent)
set(FETCHCONTENT_QUIET OFF)

set(LIBTMPL_BUILD_EXAMPLES OFF)
set(LIBTMPL_BUILD_TESTS OFF)
set(LIBTMPL_USE_ASM ON)
set(LIBTMPL_USE_MATH ON)
set(LIBTMPL_INSTALL OFF)
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

FetchContent_Declare(
    libtmpl
    GIT_REPOSITORY https://github.com/tcumby/libtmpl.git
    GIT_TAG add-cmake-build-system-redo
    SOURCE_DIR libtmpl
)

FetchContent_MakeAvailable(libtmpl)

cmake_path(GET libtmpl_SOURCE_DIR PARENT_PATH libtmpl_PARENT_DIR)
set(libtmpl_INCLUDE_DIRS "${libtmpl_PARENT_DIR}" "${libtmpl_SOURCE_DIR}" "${libtmpl_SOURCE_DIR}/include"
                         "${libtmpl_BINARY_DIR}/include"
)
