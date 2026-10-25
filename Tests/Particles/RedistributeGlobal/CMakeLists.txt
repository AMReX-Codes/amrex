foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources     main.cpp)
    if (NOT AMReX_GPU_BACKEND STREQUAL NONE)
      set(_input_files inputs.rt.cuda)
    else ()
      set(_input_files inputs.rt)
    endif ()

    # Fewer 3D steps keep the ctest short.
    if (D EQUAL 3)
        set(_ctest_params redistribute_global.nsteps=100)
    endif ()

    setup_test(${D} _sources _input_files
       CMDLINE_PARAMS_AFTER_INPUTS ${_ctest_params})

    unset(_sources)
    unset(_ctest_params)
    unset(_input_files)
endforeach()
