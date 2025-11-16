foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 2)
        return()
    endif()

    if (D EQUAL 3)
        return()
    endif()

    set(_sources  main.cpp)

    set(_input_files)

    setup_test(${D} _sources _input_files)

    unset(_sources)
    unset(_input_files)
endforeach()
