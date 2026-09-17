foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources compatibility.cpp)
    set(_input_files)
    setup_test(${D} _sources _input_files BASE_NAME Amr_GridCreation_compatibility)
    if (AMReX_AMRLEVEL)
        target_compile_definitions(Test_Amr_GridCreation_compatibility_${D}d PRIVATE TEST_AMRLEVEL)
    endif ()

    if (D EQUAL 1)
       continue()
    endif ()

    set(_sources main.cpp)
    set(_input_files inputs inputs-bf1 inputs-erf inputs-rr2 inputs-rr3-bf8 inputs-rr3-periodic
                     inputs-nochop inputs-nochop-whole inputs-nochop-rr3 inputs-nochop-periodic inputs-nochop-thin inputs-nochop-mgs0
                     inputs-rr2-bf32 inputs-mgs48 inputs-erf-nochop)

    setup_test(${D} _sources _input_files)

    # More inputs files as separate tests using the same executable.
    set(_extra inputs-rr2-bf32 inputs-rr3-bf8 inputs-mgs48)
    if (D EQUAL 3)
        list(APPEND _extra inputs-nochop inputs-nochop-whole inputs-nochop-thin
                           inputs-nochop-periodic inputs-nochop-mgs0)
    endif ()
    set(_exe ${CMAKE_CURRENT_BINARY_DIR}/${D}d/Test_Amr_GridCreation_${D}d)
    foreach(_inp IN LISTS _extra)
        set(_tname Amr_GridCreation_${_inp}_${D}d)
        if (AMReX_OMP)
            add_test(NAME ${_tname} COMMAND ${_exe} ${_inp}
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${D}d)
            set_tests_properties(${_tname} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=2)
        elseif (AMReX_MPI)
            add_test(NAME ${_tname} COMMAND mpiexec -n 2 ${_exe} ${_inp}
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${D}d)
        else ()
            add_test(NAME ${_tname} COMMAND ${_exe} ${_inp}
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${D}d)
        endif ()
    endforeach()

    unset(_sources)
    unset(_input_files)
    unset(_extra)
    unset(_exe)
endforeach()
