foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
       continue()
    endif ()

    set(_sources
       main.cpp
       MyTest.cpp
       initProb.cpp
       MyTestPlotfile.cpp
       MyTest.H
       initProb_K.H)

    set(_input_files  inputs-rt-poisson-lev )

    # A smaller 3D problem keeps the ctest short.
    if (D EQUAL 3)
        set(_ctest_params n_cell=48 max_grid_size=24)
        set(_mixed_max_error 2.5e-3)
        set(_robin_max_error 4.e-3)
        set(_com_max_error 5.5e-3)
    else ()
        set(_mixed_max_error 5.e-4)
        set(_robin_max_error 7.e-4)
        set(_com_max_error 8.e-4)
    endif ()

    setup_test(${D} _sources _input_files
       CMDLINE_PARAMS_AFTER_INPUTS ${_ctest_params})

    # The same multi-level Poisson problem with single precision multigrid
    # cycles; the error must match the double precision run.
    setup_test(${D} _sources _input_files
       BASE_NAME LinearSolvers_ABecLaplacian_C_mixed
       RUNTIME_SUBDIR mixed
       CMDLINE_PARAMS_AFTER_INPUTS ${_ctest_params} multigrid_types=geometric
                                   use_algmg_bottom=0 mixed_precision=1 max_error=${_mixed_max_error})

    # Robin BC with single precision multigrid cycles.
    set(_input_files  inputs-robin )

    setup_test(${D} _sources _input_files
       BASE_NAME LinearSolvers_ABecLaplacian_C_robin
       RUNTIME_SUBDIR robin
       CMDLINE_PARAMS_AFTER_INPUTS ${_ctest_params} mixed_precision=1 max_error=${_robin_max_error})

    # A composite solve with variable coefficients and single precision
    # multigrid cycles.  The error bound is for double precision.
    if (NOT AMReX_PRECISION STREQUAL "SINGLE")
        set(_input_files  inputs-rt-abeclap-com )

        setup_test(${D} _sources _input_files
           BASE_NAME LinearSolvers_ABecLaplacian_C_com_mixed
           RUNTIME_SUBDIR com_mixed
           CMDLINE_PARAMS_AFTER_INPUTS ${_ctest_params} multigrid_types=geometric
                                       use_algmg_bottom=0 mixed_precision=1 max_error=${_com_max_error})
    endif ()

    # Semi-coarsening exercises the line-solve smoother, which is CPU only.
    set(_input_files  inputs-rt-abeclap-semicoarsening )

    if (AMReX_GPU_BACKEND STREQUAL NONE)
        set(_no_ctest)
    else ()
        set(_no_ctest NO_CTEST)
    endif ()

    setup_test(${D} _sources _input_files
       BASE_NAME LinearSolvers_ABecLaplacian_C_semicoarsening
       RUNTIME_SUBDIR semicoarsening ${_no_ctest})

    unset(_no_ctest)
    unset(_sources)
    unset(_ctest_params)
    unset(_mixed_max_error)
    unset(_robin_max_error)
    unset(_com_max_error)
    unset(_input_files)
endforeach()
