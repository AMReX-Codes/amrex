foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources main.cpp)
    set(_input_files)

    # Fewer random domains keep the ctest short.
    setup_test(${D} _sources _input_files
       CMDLINE_PARAMS_AFTER_INPUTS ndomains=2000)

    unset(_sources)
    unset(_input_files)
endforeach()
