#
# List of subdirectories to search for CMakeLists.
#
# CUrrently missing C_BaseLib, GPU, IOBENCHMARK
set( AMREX_TESTS_SUBDIRS
   AsyncOut BaseFabTesting BBIOBenchmark complementIn FillBoundaryComparison
    )

# if (ENABLE_ASCENT)
#    list(APPEND AMREX_TESTS_SUBDIRS Blueprint)
# endif ()

# if (ENABLE_EB)
#    list(APPEND AMREX_TESTS_SUBDIRS Algoim)
# endif ()

# if (ENABLE_LINEAR_SOLVERS)
#    list(APPEND AMREX_TESTS_SUBDIRS LinearSolvers)
# endif ()

# if (ENABLE_PARTICLES)
#    list(APPEND AMREX_TESTS_SUBDIRS Particles HDF5Benchmark)
# endif ()

# if (ENABLE_SUNDIALS)
#    list(APPEND AMREX_TESTS_SUBDIRS CVODE SUNDIALS)
# endif ()

# if (ENABLE_FORTRAN_INTERFACES)
#    list(APPEND AMREX_TESTS_SUBDIRS FortranInterface)
# endif ()

# if (ENABLE_CUDA)
#    list(APPEND AMREX_TESTS_SUBDIRS GPU)
# endif ()

if (ENABLE_AMRDATA)
   list(APPEND AMREX_TESTS_SUBDIRS DataServicesTest0)
endif ()

list(TRANSFORM AMREX_TESTS_SUBDIRS PREPEND "${CMAKE_CURRENT_LIST_DIR}/")

#
# Function to setup the tutorials
#
function (setup_test _srcs  _inputs)

   cmake_parse_arguments( "" "HAS_FORTRAN_MODULES"
      "BASE_NAME;RUNTIME_SUBDIR;EXTRA_DEFINITIONS;CMDLINE_PARAMS" "" ${ARGN} )

   if (_BASE_NAME)
      set(_base_name ${_BASE_NAME})
   else ()
      string(REGEX REPLACE ".*Tests/" "" _base_name ${CMAKE_CURRENT_LIST_DIR})
      string(REPLACE "/" "_" _base_name ${_base_name})
   endif()

   if (_RUNTIME_SUBDIR)
      set(_exe_dir ${CMAKE_CURRENT_BINARY_DIR}/${_RUNTIME_SUBDIR})
   else ()
      set(_exe_dir ${CMAKE_CURRENT_BINARY_DIR})
   endif ()

   set( _exe_name  ${_base_name} )
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
         ${CMAKE_CURRENT_BINARY_DIR}/mod_files)
      set_target_properties( ${_exe_name}
         PROPERTIES
         Fortran_MODULE_DIRECTORY
         ${CMAKE_CURRENT_BINARY_DIR}/mod_files )
   endif ()

   target_link_libraries( ${_exe_name} amrex )

   if (ENABLE_CUDA)
      setup_target_for_cuda_compilation( ${_exe_name} )
   endif ()

   #
   # Assemble the commands sequence to launch the test
   #
   set(_cmd ${_test_name})

   if (_CMDLINE_PARAMS)
      list(APPEND _cmd ${_CMDLINE_PARAMS})
   endif ()

   if (${_inputs})
      file( COPY ${${_inputs}} DESTINATION ${_exe_dir} )
      list(APPEND _cmd ${${_inputs}})
   endif ()

   #
   # Add the test
   #
   add_test(
      NAME               ${_test_name}
      COMMAND            ${_cmd}
      WORKING_DIRECTORY  ${_exe_dir}
      )

endfunction ()


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
