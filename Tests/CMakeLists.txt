#
# Function to setup the tutorials
#
function (setup_test _dim _srcs  _inputs)

   cmake_parse_arguments( "" "HAS_FORTRAN_MODULES"
      "BASE_NAME;RUNTIME_SUBDIR;EXTRA_DEFINITIONS;CMDLINE_PARAMS;NTASKS;NTHREADS" "" ${ARGN} )

   if (_BASE_NAME)
      set(_base_name ${_BASE_NAME}_${_dim}d)
   else ()
      string(REGEX REPLACE ".*Tests/" "" _base_name ${CMAKE_CURRENT_LIST_DIR})
      string(REPLACE "/" "_" _base_name ${_base_name})
      set(_base_name ${_base_name}_${_dim}d)
   endif()

   if (_RUNTIME_SUBDIR)
      set(_exe_dir ${CMAKE_CURRENT_BINARY_DIR}/${_RUNTIME_SUBDIR}_${_dim}d)
   else ()
      set(_exe_dir ${CMAKE_CURRENT_BINARY_DIR}/${_dim}d)
   endif ()

   set( _exe_name  Test_${_base_name} )
   set( _test_name ${_base_name} )

   add_executable( ${_exe_name} )
   target_sources( ${_exe_name} PRIVATE ${${_srcs}} )
   set_target_properties( ${_exe_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${_exe_dir} )

   if (_EXTRA_DEFINITIONS)
      target_compile_definitions(${_exe_name} PRIVATE ${_EXTRA_DEFINITIONS})
   endif ()

   # Find out which include directory is needed
   set(_includes ${${_srcs}})
   list(FILTER _includes INCLUDE REGEX "\\.H")
   foreach(_item IN LISTS _includes)
      get_filename_component( _include_dir ${_item} DIRECTORY )
      target_include_directories( ${_exe_name} PRIVATE  ${_include_dir} )
   endforeach()

   if (_HAS_FORTRAN_MODULES)
      target_include_directories(${_exe_name}
         PRIVATE
         ${CMAKE_CURRENT_BINARY_DIR}/${_exe_name}_mod_files_${_dim}d)
      set_target_properties( ${_exe_name}
         PROPERTIES
         Fortran_MODULE_DIRECTORY
         ${CMAKE_CURRENT_BINARY_DIR}/${_exe_name}_mod_files_${_dim}d)
   endif ()

   target_link_libraries( ${_exe_name} AMReX::amrex_${_dim}d )

   if (AMReX_CUDA)
      setup_target_for_cuda_compilation( ${_exe_name} )
   endif ()

   #
   # Assemble the commands sequence to launch the test
   #
   set(_cmd ${_exe_dir}/${_exe_name})

   if (_CMDLINE_PARAMS)
      list(APPEND _cmd ${_CMDLINE_PARAMS})
   endif ()

   if (${_inputs})
      file( COPY ${${_inputs}} DESTINATION ${_exe_dir} )
      list(GET ${_inputs} 0 _first_inputs)
      get_filename_component( _inputs_filename ${_first_inputs} NAME )
      list(APPEND _cmd ${_inputs_filename})
   endif ()

   #
   # Add the test
   #
   if (AMReX_OMP)
      add_test(
         NAME               ${_test_name}
         COMMAND            ${_cmd}
         WORKING_DIRECTORY  ${_exe_dir}
      )
      set_tests_properties(${_test_name} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=2)
   elseif (AMReX_MPI)
      add_test(
         NAME               ${_test_name}
         COMMAND            mpiexec -n 2 ${_cmd}
         WORKING_DIRECTORY  ${_exe_dir}
      )
   else ()
      add_test(
         NAME               ${_test_name}
         COMMAND            ${_cmd}
         WORKING_DIRECTORY  ${_exe_dir}
      )
   endif ()

endfunction ()

if (AMReX_TEST_TYPE STREQUAL "Small")

   add_subdirectory("Amr/Advection_AmrCore")

   if (AMReX_OMP)
      add_subdirectory("OpenMP/atomicAdd")
   endif ()

   if (AMReX_PARTICLES)
      add_subdirectory("Particles/Redistribute")
   endif ()

   if (AMReX_EB)
      add_subdirectory("EB/CNS")
   endif()

   if (AMReX_LINEAR_SOLVERS)
      add_subdirectory("LinearSolvers/ABecLaplacian_C")
   endif()

else()
   #
   # List of subdirectories to search for CMakeLists.
   #
   set( AMREX_TESTS_SUBDIRS Amr AsyncOut CLZ CTOParFor DeviceGlobal Enum
                            MultiBlock MultiPeriod Parser Parser2 Reinit
                            RoundoffDomain)

   if (AMReX_PARTICLES)
     list(APPEND AMREX_TESTS_SUBDIRS Particles)
   endif ()

   if (AMReX_EB)
      list(APPEND AMREX_TESTS_SUBDIRS EB)
   endif ()

   if (AMReX_LINEAR_SOLVERS)
      list(APPEND AMREX_TESTS_SUBDIRS LinearSolvers)
   endif ()

   if (AMReX_HDF5)
      list(APPEND AMREX_TESTS_SUBDIRS HDF5Benchmark)
   endif ()

   if (AMReX_FORTRAN_INTERFACES)
      list(APPEND AMREX_TESTS_SUBDIRS FortranInterface)
   endif ()

   if (NOT AMReX_GPU_BACKEND STREQUAL NONE)
      list(APPEND AMREX_TESTS_SUBDIRS GPU)
   endif ()

   if (AMReX_GPU_BACKEND STREQUAL NONE)
      list(APPEND AMREX_TESTS_SUBDIRS OpenMP)
   endif ()

   list(TRANSFORM AMREX_TESTS_SUBDIRS PREPEND "${CMAKE_CURRENT_LIST_DIR}/")

   #
   # Loop over subdirs and add to the build those containing CMakeLists.txt
   #
   foreach (_subdir IN LISTS AMREX_TESTS_SUBDIRS)

      file( GLOB_RECURSE _tests "${_subdir}/*CMakeLists.txt" )

      foreach ( _item  IN LISTS _tests)
         get_filename_component(_dir ${_item} DIRECTORY )
         add_subdirectory(${_dir})
      endforeach ()

   endforeach ()
endif()
