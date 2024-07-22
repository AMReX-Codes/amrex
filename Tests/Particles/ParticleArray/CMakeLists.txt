foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources main.cpp)

    setup_test(${D} _sources FALSE)

    unset(_sources)
endforeach()
