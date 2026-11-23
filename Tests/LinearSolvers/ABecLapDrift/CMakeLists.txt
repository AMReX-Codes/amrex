# The tolerances are for double precision
if (AMReX_PRECISION STREQUAL "SINGLE")
    return()
endif ()

foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
       continue()
    endif ()

    set(_sources main.cpp)
    set(_input_files inputs)

    if (D EQUAL 3)
        set(_ctest_params n_cell=24 max_grid_size=12)
    endif ()

    # One level: apply against the operator written out, and the solve
    setup_test(${D} _sources _input_files
       CMDLINE_PARAMS_AFTER_INPUTS ${_ctest_params})

    # Composite solve over two levels, with Jacobi smoothing
    setup_test(${D} _sources _input_files
       BASE_NAME LinearSolvers_ABecLapDrift_amr
       RUNTIME_SUBDIR amr
       CMDLINE_PARAMS_AFTER_INPUTS ${_ctest_params} max_level=1 use_gauss_seidel=0)

    # Without drift: MLABecLaplacian
    setup_test(${D} _sources _input_files
       BASE_NAME LinearSolvers_ABecLapDrift_nodrift
       RUNTIME_SUBDIR nodrift
       CMDLINE_PARAMS_AFTER_INPUTS ${_ctest_params} cscalar=0)

    unset(_sources)
    unset(_input_files)
    unset(_ctest_params)
endforeach()
