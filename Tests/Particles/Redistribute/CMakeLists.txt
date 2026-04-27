foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources     main.cpp)
    if (NOT AMReX_GPU_BACKEND STREQUAL NONE)
      set(_input_files inputs.rt.cuda  )
    else ()
      set(_input_files inputs.rt  )
    endif ()

    setup_test(${D} _sources _input_files)

    if (D EQUAL 3 AND NOT AMReX_GPU_BACKEND STREQUAL NONE)
      set(_exe_dir ${CMAKE_CURRENT_BINARY_DIR}/${D}d)
      set(_exe_name Test_Particles_Redistribute_${D}d)

      file(COPY inputs.rt.cuda.tiling DESTINATION ${_exe_dir})

      if (AMReX_OMP)
        add_test(
          NAME              Particles_Redistribute_${D}d_Tiling
          COMMAND           ${_exe_dir}/${_exe_name} inputs.rt.cuda.tiling
          WORKING_DIRECTORY ${_exe_dir}
        )
        set_tests_properties(Particles_Redistribute_${D}d_Tiling PROPERTIES ENVIRONMENT OMP_NUM_THREADS=2)
      elseif (AMReX_MPI)
        add_test(
          NAME              Particles_Redistribute_${D}d_Tiling
          COMMAND           mpiexec -n 2 ${_exe_dir}/${_exe_name} inputs.rt.cuda.tiling
          WORKING_DIRECTORY ${_exe_dir}
        )
      else ()
        add_test(
          NAME              Particles_Redistribute_${D}d_Tiling
          COMMAND           ${_exe_dir}/${_exe_name} inputs.rt.cuda.tiling
          WORKING_DIRECTORY ${_exe_dir}
        )
      endif ()
    endif ()

    unset(_sources)
    unset(_input_files)
endforeach()
