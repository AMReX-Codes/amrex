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

    setup_test(${D} _sources _input_files)

    # Semi-coarsening exercises the line-solve smoother, which is CPU only.
    set(_input_files  inputs-rt-abeclap-semicoarsening )

    if (AMReX_GPU_BACKEND STREQUAL NONE)
        set(_no_ctest)
    else ()
        set(_no_ctest NO_CTEST)
    endif ()

    setup_test(${D} _sources _input_files
       BASE_NAME LinearSolvers_ABecLaplacian_C_semicoarsening ${_no_ctest})

    unset(_no_ctest)
    unset(_sources)
    unset(_input_files)
endforeach()
