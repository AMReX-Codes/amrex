# MLNodeLaplacian is built with AMReX_LINEAR_SOLVERS_INCFLO
if (NOT AMReX_LINEAR_SOLVERS_INCFLO)
    return()
endif ()

foreach(D IN LISTS AMReX_SPACEDIM)
    # MLNodeLaplacian::compRHS does not support 1D
    if (D EQUAL 1)
       continue()
    endif ()

    set(_sources main.cpp checks.cpp study.cpp NodalProjectionWall.H)

    set(_input_files inputs-ci)

    setup_test(${D} _sources _input_files)

    # The same checks without the correction: the wall nodes are then first order
    setup_test(${D} _sources _input_files
       BASE_NAME LinearSolvers_NodalProjectionWall_NoCorrection
       CMDLINE_PARAMS_AFTER_INPUTS wall_corr=0)

    unset(_sources)
    unset(_input_files)
endforeach()
