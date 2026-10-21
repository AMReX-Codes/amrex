foreach(D IN LISTS AMReX_SPACEDIM)
    if ( (D EQUAL 1) OR (NOT AMReX_FORTRAN_INTERFACES) )
        continue()
    endif ()

    set(_sources main.F90 mytest.F90 init_prob.F90)

    file( GLOB_RECURSE _input_files LIST_DIRECTORIES false ${CMAKE_CURRENT_LIST_DIR}/input*)

    # A smaller 3D problem keeps the ctest short.
    if (D EQUAL 3)
        set(_ctest_params n_cell=64 max_grid_size=32)
    endif ()

    setup_test(${D} _sources _input_files HAS_FORTRAN_MODULES
       CMDLINE_PARAMS_AFTER_INPUTS ${_ctest_params})

    unset(_sources)
    unset(_ctest_params)
    unset(_input_files)
endforeach()
