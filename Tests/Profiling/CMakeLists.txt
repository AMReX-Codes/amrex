foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources main.cpp)
    set(_input_files inputs)

    setup_test(${D} _sources _input_files)

    unset(_sources)
    unset(_input_files)

    # Same executable with all timers synchronized: the report must stay complete.
    set(_cmd ${CMAKE_CURRENT_BINARY_DIR}/${D}d/Test_Profiling_${D}d inputs
        tiny_profiler.device_synchronize_around_region=1
        tiny_profiler.output_file=tiny_profiler_device_sync.out)
    if (AMReX_MPI AND NOT AMReX_OMP)
        list(PREPEND _cmd mpiexec -n 2)
    endif ()
    add_test(NAME Profiling_DeviceSync_${D}d
             COMMAND ${_cmd}
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${D}d)
    if (AMReX_OMP)
        set_tests_properties(Profiling_DeviceSync_${D}d PROPERTIES ENVIRONMENT OMP_NUM_THREADS=2)
    endif ()
    unset(_cmd)
endforeach()
