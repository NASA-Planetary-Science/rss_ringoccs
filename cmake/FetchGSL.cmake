cmake_minimum_required(VERSION 3.20)
include_guard()

if(APPLE OR UNIX)
    include(FindGSL)
else()
    # vcpkg installs GSL for Windows via the vcpkg.json manifest
endif()

find_package(GSL REQUIRED)
if(NOT GSL_FOUND)
    message(FATAL_ERROR "GSL was not found")
endif()
