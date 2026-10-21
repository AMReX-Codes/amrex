if (NOT AMReX_LINEAR_SOLVERS_EM)
    return()
endif ()

foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources
       main.cpp
       MyTest.cpp
       MyTest.H
       initProb.cpp
       initProb_K.H)

    set(_input_files  inputs)

    # A smaller 3D problem keeps the ctest short.
    if (D EQUAL 3)
        set(_ctest_params n_cell=32 max_grid_size=16)
    endif ()

    setup_test(${D} _sources _input_files
       CMDLINE_PARAMS_AFTER_INPUTS ${_ctest_params})

    unset(_sources)
    unset(_ctest_params)
    unset(_input_files)
endforeach()
