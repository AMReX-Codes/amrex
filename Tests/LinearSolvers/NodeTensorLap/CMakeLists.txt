if (NOT AMReX_LINEAR_SOLVERS_EM)
    return()
endif ()

if (AMReX_GPU_BACKEND STREQUAL NONE)
    foreach(D IN LISTS AMReX_SPACEDIM)
        if (D EQUAL 1)
            continue()
        endif ()

        set(_sources main.cpp MyTest.cpp MyTest.H MyTestPlotfile.cpp)
        set(_input_files)

        # A smaller 3D problem keeps the ctest short.
        if (D EQUAL 3)
            set(_ctest_params n_cell=64 max_grid_size=32)
        endif ()

        setup_test(${D} _sources _input_files
           CMDLINE_PARAMS "multigrid_types=algebraic hybrid geometric"
           CMDLINE_PARAMS_AFTER_INPUTS ${_ctest_params})

        unset(_sources)
        unset(_ctest_params)
        unset(_input_files)
    endforeach()
endif ()
