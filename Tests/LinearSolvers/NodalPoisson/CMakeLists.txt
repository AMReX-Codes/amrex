if (NOT AMReX_LINEAR_SOLVERS_INCFLO)
    return()
endif ()

foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
       continue()
    endif ()

    set(_sources main.cpp MyTest.cpp MyTest.H MyTestPlotfile.cpp)

    set(_input_files inputs-ci)

    setup_test(${D} _sources _input_files)

    set(_input_files inputs-sigma-update)

    setup_test(${D} _sources _input_files
       BASE_NAME LinearSolvers_NodalPoisson_SigmaUpdate)

    unset(_sources)
    unset(_input_files)
endforeach()
