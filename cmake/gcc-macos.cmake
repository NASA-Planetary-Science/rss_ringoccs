include_guard()

# Find gcc as installed by Homebrew
find_program(
    GCC
    NAMES gcc gcc-9 gcc-10 gcc-11
    PATHS /usr/local/Cellar /usr/local/opt /opt/local/bin /usr/local/bin
    NO_DEFAULT_PATH
)

if(GCC)
    set(CMAKE_C_COMPILER ${GCC})

    set(extraOpts
        -ansi
        -pedantic
        -pedantic-errors
        -Wall
        -Wextra
        -Wpedantic
        -Wmisleading-indentation
        -Wmissing-field-initializers
        -Wmissing-prototypes
        -Wold-style-definition
        -Winit-self
        -Wmissing-declarations
        -Wnull-dereference
        -Wwrite-strings
        -Wstrict-prototypes
        -fuse-linker-plugin
        -g
        -fPIC
        -O3
    )
    set(CMAKE_C_FLAGS_DEBUG_INIT ${extraOpts})

    if(CMAKE_C_COMPILE_OPTIONS_IPO MATCHES "-fno-fat-lto-objects")
        string(
            REPLACE
            "-fno-fat-lto-objects"
            ""
            ${CMAKE_C_COMPILE_OPTIONS_IPO}
            CMAKE_C_COMPILE_OPTIONS_IPO
        )
    endif()
else()
    message(FATAL_ERROR "gcc was not found")
endif()
