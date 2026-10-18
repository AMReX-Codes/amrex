if (NOT AMReX_EB)
   return()
endif ()

foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
       continue()
    endif ()

    set(_sources main.cpp)
    set(_input_files inputs)

    # A benchmark, not a regression test.
    setup_test(${D} _sources _input_files NO_CTEST)

    unset(_sources)
    unset(_input_files)
endforeach()
