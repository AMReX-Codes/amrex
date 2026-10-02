foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources main.cpp)
    set(_input_files inputs)

    if (AMReX_GPU_BACKEND STREQUAL NONE)
        setup_test(${D} _sources _input_files)
    else ()
        # Every case syncs the GPU often, which is slow when the GPU is
        # shared, so run the variations for one problem and nothing else
        # alongside.
        setup_test(${D} _sources _input_files
            CMDLINE_PARAMS_AFTER_INPUTS variation_problems=checker)
        set_tests_properties(Algebra_AlgMG_${D}d PROPERTIES RUN_SERIAL TRUE)
    endif ()

    unset(_sources)
    unset(_input_files)
endforeach()
