foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources     main.cpp)
    set(_input_files config.toml config-2.toml)

    setup_test(${D} _sources _input_files)

    unset(_sources)
    unset(_input_files)
endforeach()
