if (AMReX_GPU_BACKEND STREQUAL "NONE")
    return()
endif ()

foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources main.cpp)
    set(_input_files)

    if (AMReX_LINEAR_SOLVERS)
        setup_test(${D} _sources _input_files
                   EXTRA_DEFINITIONS AMREX_TEST_LINEAR_SOLVERS)
    else ()
        setup_test(${D} _sources _input_files)
    endif ()

    unset(_sources)
    unset(_input_files)
endforeach()
