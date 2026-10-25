foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
       continue()
    endif ()

    set(_sources  main.cpp)
    set(_input_files  inputs plot_eb.py)

    # Fewer 3D cases keep the ctest short.
    if (D EQUAL 3)
        set(_ctest_params ncases=10)
    endif ()

    setup_test(${D} _sources _input_files
       CMDLINE_PARAMS_AFTER_INPUTS ${_ctest_params})

    unset(_sources)
    unset(_ctest_params)
    unset(_input_files)
endforeach()
