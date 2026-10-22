foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources     main.cpp     mypc.H     trilinear_deposition_K.H)
    set(_input_files inputs  )

    # A smaller 3D problem keeps the ctest short.
    if (D EQUAL 3)
        set(_ctest_params nx=64 ny=64 nz=64)
    endif ()

    setup_test(${D} _sources _input_files
       CMDLINE_PARAMS_AFTER_INPUTS ${_ctest_params})

    unset(_sources)
    unset(_ctest_params)
    unset(_input_files)
endforeach()
