if (NOT AMReX_LINEAR_SOLVERS_INCFLO)
    return()
endif ()

foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
       continue()
    endif ()

    set(_sources main.cpp MyTest.cpp MyTest.H MyTestPlotfile.cpp)

    set(_input_files inputs-ci)

    # A smaller 3D problem keeps the ctest short.
    if (D EQUAL 3)
        set(_ctest_params n_cell=24 max_grid_size=12)
    endif ()

    setup_test(${D} _sources _input_files
       CMDLINE_PARAMS_AFTER_INPUTS ${_ctest_params})

    set(_input_files inputs-sigma-update)

    setup_test(${D} _sources _input_files
       BASE_NAME LinearSolvers_NodalPoisson_SigmaUpdate)

    unset(_sources)
    unset(_ctest_params)
    unset(_input_files)
endforeach()
