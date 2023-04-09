include_guard()

find_program(GCC gcc)
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
        -g
        -fPIC
        -O3
    )
    set(CMAKE_C_FLAGS_DEBUG_INIT ${extraOpts})
else()
    message(FATAL_ERROR "gcc was not found")
endif()
