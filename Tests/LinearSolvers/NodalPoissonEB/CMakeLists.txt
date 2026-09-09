if (NOT AMReX_EB OR NOT AMReX_LINEAR_SOLVERS_EM)
    return()
endif ()

foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
        continue()
    endif ()

    set(_sources main.cpp MyTest.cpp MyTest.H MyTestPlotfile.cpp)

    set(_input_files inputs-${D}d)

    setup_test(${D} _sources _input_files)

    if (D EQUAL 3)
        set(_input_files inputs-3d-stl sphere.stl)

        # The STL test does not work in single precision, so only build it there.
        if (AMReX_PRECISION STREQUAL "SINGLE")
            set(_no_ctest NO_CTEST)
        else ()
            set(_no_ctest)
        endif ()

        setup_test(${D} _sources _input_files
            BASE_NAME LinearSolvers_NodalPoissonEB_STL ${_no_ctest})

        unset(_no_ctest)
    endif ()

    unset(_sources)
    unset(_input_files)
endforeach()
