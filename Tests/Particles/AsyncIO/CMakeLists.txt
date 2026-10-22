foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources     main.cpp)
    set(_input_files inputs  )

    # Fewer 3D steps keep the ctest short.
    if (D EQUAL 3)
        set(_ctest_params nsteps=400)
    endif ()

    setup_test(${D} _sources _input_files
       CMDLINE_PARAMS_AFTER_INPUTS ${_ctest_params})

    unset(_sources)
    unset(_ctest_params)
    unset(_input_files)
endforeach()
