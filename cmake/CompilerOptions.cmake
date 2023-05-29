include_guard()

if(APPLE)
    include(gcc-macos)
elseif(UNIX)
    include(gcc)
elseif(WIN32)
    if(MSVC)
        # Compile as C++ for Microsoft Visual Studio to avoid the non-standard complex number implementation in the C
        # compiler
        list(APPEND CMAKE_CXX_SOURCE_FILE_EXTENSIONS c)
        set(CMAKE_CXX_STANDARD 11)
        set(CMAKE_GENERATOR_PLATFORM
            x64
            CACHE INTERNAL ""
        )

        # Enable parallel builds
        add_compile_options("/MP")
    endif()
endif()
