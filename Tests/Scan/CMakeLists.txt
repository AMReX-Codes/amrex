foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources main.cpp)
    set(_input_files)

    setup_test(${D} _sources _input_files)
    if (AMReX_SYCL)
        setup_test(${D} _sources _input_files
            BASE_NAME ScanSinglePass
            RUNTIME_SUBDIR SinglePass
            EXTRA_DEFINITIONS AMREX_SYCL_NO_MULTIPASS_SCAN)
    endif ()

    unset(_sources)
    unset(_input_files)
endforeach()
