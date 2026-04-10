foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources main.cpp)

    setup_test(${D} _sources _input_files
        BASE_NAME GPU_ArenaLifetime)

    unset(_sources)
    unset(_input_files)
endforeach()
