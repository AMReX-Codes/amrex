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
    set(_input_files inputs inputs-bf1 inputs-rr2 inputs-rr3-bf8 inputs-rr3-periodic
                     inputs-no-box-split inputs-no-box-split-whole inputs-no-box-split-rr3
                     inputs-no-box-split-periodic inputs-no-box-split-thin inputs-no-box-split-mgs0
                     inputs-rr2-bf32 inputs-mgs48 inputs-erf-no-box-split
                     inputs-no-box-split-merge-thin inputs-no-box-split-mgs-unit inputs-no-box-split-nesting)

    # In 2D, inputs files with no_box_split_dir = 2 use direction 1 instead, and
    # 2D tag centers so that tagged regions do not share columns.
    set(_args)
    if (D EQUAL 2)
        set(_args amr.no_box_split_dir=1 test.centers=0.3 0.4 0.9 0.95)
    endif ()

    setup_test(${D} _sources _input_files CMDLINE_PARAMS_AFTER_INPUTS ${_args})

    # More inputs files as separate tests using the same executable.
    set(_extra inputs-rr2-bf32 inputs-rr3-bf8 inputs-mgs48 inputs-no-box-split-mgs-unit)
    if (D EQUAL 3)
        list(APPEND _extra inputs-no-box-split inputs-no-box-split-whole inputs-no-box-split-thin
                           inputs-no-box-split-periodic inputs-no-box-split-mgs0
                           inputs-no-box-split-merge-thin inputs-no-box-split-nesting)
    endif ()
    set(_exe ${CMAKE_CURRENT_BINARY_DIR}/${D}d/Test_Amr_GridCreation_${D}d)
    foreach(_inp IN LISTS _extra)
        set(_inp_args)
        file(STRINGS ${CMAKE_CURRENT_LIST_DIR}/${_inp} _no_box_split REGEX "no_box_split_dir")
        if (_no_box_split)
            set(_inp_args ${_args})
        endif ()
        set(_tname Amr_GridCreation_${_inp}_${D}d)
        if (AMReX_OMP)
            add_test(NAME ${_tname} COMMAND ${_exe} ${_inp} ${_inp_args}
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${D}d)
            set_tests_properties(${_tname} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=2)
        elseif (AMReX_MPI)
            add_test(NAME ${_tname} COMMAND mpiexec -n 2 ${_exe} ${_inp} ${_inp_args}
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${D}d)
        else ()
            add_test(NAME ${_tname} COMMAND ${_exe} ${_inp} ${_inp_args}
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${D}d)
        endif ()
    endforeach()

    unset(_sources)
    unset(_input_files)
    unset(_extra)
    unset(_args)
    unset(_inp_args)
    unset(_no_box_split)
    unset(_exe)
endforeach()
