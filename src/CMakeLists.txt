cmake_minimum_required(VERSION 3.20)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    if(APPLE)
        include(gcc-macos)
    elseif(UNIX)
        include(gcc)
    endif()
endif()

include(GenerateExportHeader)

project(rssringoccs_library)

# if(MSVC) file(GLOB_RECURSE CFILES "*.c") set_source_files_properties(${CFILES} PROPERTIES LANGUAGE CXX) endif()
find_package(cspice REQUIRED)
add_library(librssringoccs SHARED)
if(UNIX)
    target_link_libraries(librssringoccs PRIVATE m)
endif()
find_package(OpenMP)
if(OpenMP_C_FOUND)
    message(STATUS "OpenMP ${OpenMP_C_VERSION} found")
    target_link_libraries(librssringoccs PUBLIC OpenMP::OpenMP_C)
endif()

if(IPO_SUPPORTED)
    set_property(TARGET librssringoccs PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(WARNING "IPO is not supported: ${output}")
endif()

target_include_directories(
    librssringoccs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include ${CMAKE_CURRENT_SOURCE_DIR}/../..
                          ${libtmpl_SOURCE_DIR}/..
)

add_subdirectory("csv_tools")
add_subdirectory("fresnel_transform")
add_subdirectory("history")
add_subdirectory("occultation_geometry")
add_subdirectory("reconstruction")

target_compile_definitions(librssringoccs PUBLIC -DNDEBUG)
generate_export_header(
    librssringoccs # generates the export header shared_EXPORTS.h automatically
    BASE_NAME
    librssringoccs
    EXPORT_MACRO_NAME
    RSS_RINGOCCS_EXPORT
    EXPORT_FILE_NAME
    librssringoccs_exports.h
    STATIC_DEFINE
    SHARED_EXPORTS_BUILT_AS_STATIC
)
target_include_directories(librssringoccs PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(librssringoccs PRIVATE tmpl::tmpl cspice::cspice)
file(GLOB LIBRSSRINGOCCS_HEADERS ../include/*.h)
set_target_properties(
    librssringoccs
    PROPERTIES LANGUAGE $<IF:$<C_COMPILER_ID:MSVC>CXX,C>
               INTERPROCEDURAL_OPTIMIZATION ${IPO_SUPPORTED}
               PUBLIC_HEADER "${LIBRSSRINGOCCS_HEADERS}"
               OUTPUT_NAME rssringoccs
               IMPORTED_LOCATION $<TARGET_FILE:librssringoccs>
               IMPORTED_IMPLIB $<TARGET_LINKER_FILE:librssringoccs>
               INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_BINARY_DIR}
)

add_library(rss::librssringoccs ALIAS librssringoccs)

if(NOT DEFINED SKBUILD)
    include(GNUInstallDirs)
    install(
        TARGETS librssringoccs
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    )
    install(
        TARGETS tmpl
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
else()
    if(MACOSX)
        set(LIBRSSRINGOCCS_LIB_PACKAGE_DIR
            ".dylibs"
            CACHE INTERNAL "The folder in which to store \
        librssringoccs.so/.dylib/.dll"
        )
    else()
        set(LIBRSSRINGOCCS_LIB_PACKAGE_DIR
            ".libs"
            CACHE INTERNAL "The folder in which to store \
        librssringoccs.so/.dylib/.dll"
        )
    endif()

    configure_file("../cmake/__config__.py.in" "__config__.py")

    install(
        TARGETS librssringoccs tmpl
        LIBRARY DESTINATION ${LIBRSSRINGOCCS_LIB_PACKAGE_DIR}
        RUNTIME DESTINATION ${LIBRSSRINGOCCS_LIB_PACKAGE_DIR}
    )
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/__config__.py" DESTINATION .)
endif()
