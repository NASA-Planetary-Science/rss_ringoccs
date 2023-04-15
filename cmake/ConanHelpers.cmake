function(conan_install)
    set(prefix _CONAN_INSTALL)
    set(noValues "")
    set(singleValues SOURCE_DIR OUTPUT_DIR)
    set(multiValues "")
    # Process the arguments passed in
    include(CMakeParseArguments)
    cmake_parse_arguments(${prefix} "${noValues}" "${singleValues}" "${multiValues}" ${ARGN})

    execute_process(
        COMMAND conan install --build=missing --output-folder ${${prefix}OUTPUT_DIR} ${${prefix}SOURCE_DIR}
        COMMAND_ECHO STDOUT
        COMMAND_ERROR_IS_FATAL ANY
    )
endfunction()
