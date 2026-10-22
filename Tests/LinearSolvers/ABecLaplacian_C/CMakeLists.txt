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
    endif ()

    setup_test(${D} _sources _input_files
       CMDLINE_PARAMS_AFTER_INPUTS ${_ctest_params})

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
    unset(_input_files)
endforeach()
