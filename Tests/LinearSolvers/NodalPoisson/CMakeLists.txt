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
        set(_mixed_max_error 3.5e-2)
    else ()
        set(_mixed_max_error 3.e-3)
    endif ()

    setup_test(${D} _sources _input_files
       CMDLINE_PARAMS_AFTER_INPUTS ${_ctest_params})

    # Single precision multigrid cycles with an overset mask on both AMR
    # levels; the error must match the double precision run.
    setup_test(${D} _sources _input_files
       BASE_NAME LinearSolvers_NodalPoisson_mixed
       RUNTIME_SUBDIR mixed
       CMDLINE_PARAMS_AFTER_INPUTS ${_ctest_params} multigrid_types=geometric
                                   mixed_precision=1 overset=1 max_error=${_mixed_max_error})

    set(_input_files inputs-sigma-update)

    setup_test(${D} _sources _input_files
       BASE_NAME LinearSolvers_NodalPoisson_SigmaUpdate)

    setup_test(${D} _sources _input_files
       BASE_NAME LinearSolvers_NodalPoisson_SigmaUpdate_mixed
       RUNTIME_SUBDIR sigma_update_mixed
       CMDLINE_PARAMS_AFTER_INPUTS mixed_precision=1)

    unset(_sources)
    unset(_ctest_params)
    unset(_mixed_max_error)
    unset(_input_files)
endforeach()
