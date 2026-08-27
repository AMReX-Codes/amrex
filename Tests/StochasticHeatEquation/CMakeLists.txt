foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
        set(_sources main.cpp)
        set(_input_files inputs-ci)

        setup_test(${D} _sources _input_files BASE_NAME StochasticHeatEquation)

        unset(_sources)
        unset(_input_files)
    endif()
endforeach()
