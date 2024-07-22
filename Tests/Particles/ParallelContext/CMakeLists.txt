foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources     main.cpp)
    if (NOT AMReX_GPU_BACKEND STREQUAL NONE)
      set(_input_files inputs.rt.cuda  )
    else ()
      set(_input_files inputs.rt  )
    endif ()

    setup_test(${D} _sources _input_files)

    unset(_sources)
    unset(_input_files)
endforeach()
