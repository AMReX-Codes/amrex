foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources main.cpp)
    set(_input_files inputs)

    if (AMReX_GPU_BACKEND STREQUAL NONE)
        setup_test(${D} _sources _input_files)
    else ()
        # Every case syncs the GPU often, which is slow when the GPU is
        # shared, so run the variations for one problem, split them across
        # dimensions, and run nothing else alongside.
        if (D EQUAL 1)
            set(_groups bottom singular)
        elseif (D EQUAL 2)
            set(_groups krylov smoother)
        else ()
            set(_groups interp coarsening truncation)
        endif ()
        setup_test(${D} _sources _input_files
            CMDLINE_PARAMS_AFTER_INPUTS variation_problems=checker variation_groups=${_groups})
        set_tests_properties(Algebra_AlgMG_${D}d PROPERTIES RUN_SERIAL TRUE)
        unset(_groups)
    endif ()

    unset(_sources)
    unset(_input_files)
endforeach()
