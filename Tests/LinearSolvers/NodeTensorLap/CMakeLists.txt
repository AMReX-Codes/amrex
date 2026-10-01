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

        setup_test(${D} _sources _input_files
           CMDLINE_PARAMS "multigrid_types=algebraic hybrid geometric")

        unset(_sources)
        unset(_input_files)
    endforeach()
endif ()
