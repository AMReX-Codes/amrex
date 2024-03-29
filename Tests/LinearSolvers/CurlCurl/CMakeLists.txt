foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
       return()
    endif ()

    set(_sources
       main.cpp
       MyTest.cpp
       MyTest.H
       initProb.cpp
       initProb_K.H)

    set(_input_files  inputs)

    setup_test(${D} _sources _input_files)

    unset(_sources)
    unset(_input_files)
endforeach()
