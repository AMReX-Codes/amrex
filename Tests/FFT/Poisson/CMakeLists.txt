foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources  main.cpp)

    set(_input_files)

    # A smaller 3D problem keeps the ctest short.
    if (D EQUAL 3)
        set(_ctest_params n_cell_x=32 n_cell_y=24 n_cell_z=64)
    endif ()

    setup_test(${D} _sources _input_files
       CMDLINE_PARAMS_AFTER_INPUTS ${_ctest_params})

    unset(_sources)
    unset(_ctest_params)
    unset(_input_files)
endforeach()
