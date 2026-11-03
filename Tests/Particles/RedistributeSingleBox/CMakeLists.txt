foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources     main.cpp)
    set(_input_files )

    setup_test(${D} _sources _input_files)

    # setup_test runs MPI builds on two processes, the in-place Redistribute needs one
    if (AMReX_MPI AND NOT AMReX_OMP)
        set(_exe_dir ${CMAKE_CURRENT_BINARY_DIR}/${D}d)
        add_test(
            NAME              Particles_RedistributeSingleBox_${D}d_OneProcess
            COMMAND           mpiexec -n 1 ${_exe_dir}/Test_Particles_RedistributeSingleBox_${D}d
            WORKING_DIRECTORY ${_exe_dir}
        )
        unset(_exe_dir)
    endif ()

    unset(_sources)
    unset(_input_files)
endforeach()
