if (( (AMReX_GPU_BACKEND STREQUAL "CUDA") OR
      (AMReX_GPU_BACKEND STREQUAL "HIP" ) ) AND
    (NOT AMReX_GPU_RDC))
   return()
endif ()

foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources     main.cpp global_vars.cpp init.cpp work.cpp)
    set(_input_files)

    setup_test(${D} _sources _input_files)

    unset(_sources)
    unset(_input_files)
endforeach()
