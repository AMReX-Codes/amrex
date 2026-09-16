foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
       continue()
    endif ()

    set(_sources main.cpp)
    set(_input_files inputs inputs-bf1 inputs-erf inputs-rr2 inputs-rr3-bf8 inputs-rr3-periodic
                     inputs-nochop inputs-nochop-whole inputs-nochop-rr3 inputs-nochop-periodic inputs-erf-nochop)

    setup_test(${D} _sources _input_files)

    unset(_sources)
    unset(_input_files)
endforeach()
