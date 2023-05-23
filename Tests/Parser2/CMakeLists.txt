foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources     main.cpp fn.cpp)
    set(_input_files fn.cpp)

    setup_test(${D} _sources _input_files)

    unset(_sources)
    unset(_input_files)
endforeach()
