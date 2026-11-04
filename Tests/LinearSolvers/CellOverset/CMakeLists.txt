foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
       continue()
    endif ()

    set(_sources main.cpp MyTest.cpp MyTest.H)
    set(_input_files inputs)

    setup_test(${D} _sources _input_files
       CMDLINE_PARAMS_AFTER_INPUTS ${AMREX_TEST_MLMG_TYPE_PARAMS})

    unset(_sources)
    unset(_input_files)
endforeach()
