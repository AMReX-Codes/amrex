foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources     main.cpp)
    set(_input_files inputs)

    setup_test(${D} _sources _input_files
       CMDLINE_PARAMS_AFTER_INPUTS unset_cmdline_test=1 UNSET=unset_from_cmdline)

    unset(_sources)
    unset(_input_files)
endforeach()
