foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources main.cpp)
    set(_input_files inputs inputs_hide_async inputs_other)

    setup_test(${D} _sources _input_files BASE_NAME TinyProfilerAsync)

    set(_exe_dir ${CMAKE_CURRENT_BINARY_DIR}/${D}d)
    set(_exe_name Test_TinyProfilerAsync_${D}d)
    set(_hidden_test TinyProfilerAsyncHidden_${D}d)
    set(_other_test TinyProfilerAsyncOther_${D}d)

    if (AMReX_OMP)
        add_test(NAME ${_hidden_test}
                 COMMAND ${_exe_dir}/${_exe_name} inputs_hide_async
                 WORKING_DIRECTORY ${_exe_dir})
        set_tests_properties(${_hidden_test} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=2)
    elseif (AMReX_MPI)
        add_test(NAME ${_hidden_test}
                 COMMAND mpiexec -n 2 ${_exe_dir}/${_exe_name} inputs_hide_async
                 WORKING_DIRECTORY ${_exe_dir})
    else ()
        add_test(NAME ${_hidden_test}
                 COMMAND ${_exe_dir}/${_exe_name} inputs_hide_async
                 WORKING_DIRECTORY ${_exe_dir})
    endif ()

    if (AMReX_OMP)
        add_test(NAME ${_other_test}
                 COMMAND ${_exe_dir}/${_exe_name} inputs_other
                 WORKING_DIRECTORY ${_exe_dir})
        set_tests_properties(${_other_test} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=2)
    elseif (AMReX_MPI)
        add_test(NAME ${_other_test}
                 COMMAND mpiexec -n 2 ${_exe_dir}/${_exe_name} inputs_other
                 WORKING_DIRECTORY ${_exe_dir})
    else ()
        add_test(NAME ${_other_test}
                 COMMAND ${_exe_dir}/${_exe_name} inputs_other
                 WORKING_DIRECTORY ${_exe_dir})
    endif ()

    unset(_sources)
    unset(_input_files)
endforeach()
