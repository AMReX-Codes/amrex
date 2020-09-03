#
# List of subdirectories to search for CMakeLists.
# For now, we do not include BluePrint, FOrkJoin, GPU, MUI, SDC, SWFFT
#
set( AMREX_TUTORIALS_SUBDIRS Amr Basic )

if (ENABLE_EB)
   list(APPEND AMREX_TUTORIALS_SUBDIRS EB)
endif ()

if (ENABLE_LINEAR_SOLVERS)
   list(APPEND AMREX_TUTORIALS_SUBDIRS LinearSolvers)
endif ()

if (ENABLE_PARTICLES)
   list(APPEND AMREX_TUTORIALS_SUBDIRS Particles)
endif ()

if (ENABLE_SUNDIALS)
   list(APPEND AMREX_TUTORIALS_SUBDIRS CVODE SUNDIALS)
endif ()

if (ENABLE_FORTRAN_INTERFACES)
   list(APPEND AMREX_TUTORIALS_SUBDIRS FortranInterface)
endif ()

list(TRANSFORM AMREX_TUTORIALS_SUBDIRS PREPEND "${CMAKE_CURRENT_LIST_DIR}/")


#
# Macro to setup the tutorials
#
macro (setup_tutorial _sources  _input_files)

   cmake_parse_arguments( "" "HAS_FORTRAN_MODULES" "BASE_NAME;RUNTIME_SUBDIR" "" ${ARGN} )

   if (_BASE_NAME)
      set(_base_name ${_BASE_NAME})
      unset(_BASE_NAME)
   else ()
      get_filename_component( _base_name ${CMAKE_CURRENT_LIST_DIR} NAME )
   endif()

   if (_RUNTIME_SUBDIR)
      set(_exe_dir ${CMAKE_CURRENT_BINARY_DIR}/${_RUNTIME_SUBDIR})
      unset(_RUNTIME_SUBDIR)
   else ()
      set(_exe_dir ${CMAKE_CURRENT_BINARY_DIR})
   endif ()

   set( _test_name ${_base_name} )
   set( _exe_name  ${_base_name}.exe )

   add_executable( ${_exe_name} )

   target_sources( ${_exe_name} PRIVATE ${${_sources}} )

   set_target_properties( ${_exe_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${_exe_dir} )

   # Find out which include directory is needed
   set(_includes ${${_sources}})
   list(FILTER _includes INCLUDE REGEX "\\.H")
   foreach(_item IN LISTS _includes)
      get_filename_component( _include_dir ${_item} DIRECTORY)
      target_include_directories( ${_exe_name} PRIVATE  ${_include_dir} )
   endforeach()
   unset(_includes)

   if (ENABLE_CUDA)
      set_cpp_sources_to_cuda_language( ${_exe_name} )
   endif ()

   # Extra arguments
   if (_HAS_FORTRAN_MODULES)
      target_include_directories(${_exe_name}
         PRIVATE
         ${CMAKE_CURRENT_BINARY_DIR}/mod_files)
      set_target_properties( ${_exe_name}
         PROPERTIES
         Fortran_MODULE_DIRECTORY
         ${CMAKE_CURRENT_BINARY_DIR}/mod_files )
      unset(_HAS_FORTRAN_MODULES)
   endif ()

   target_link_libraries( ${_exe_name} amrex )

   if (${_input_files})
      file( COPY ${${_input_files}} DESTINATION ${_exe_dir} )
   endif ()

   unset(_base_name)
   unset(_test_name)
   unset(_exe_name)
   unset(_exe_dir)

endmacro ()

macro (add_tutorial targets_list)
endmacro ()


#
# Loop over subdirs and add to the build those containing CMakeLists.txt
#
foreach (_subdir IN LISTS AMREX_TUTORIALS_SUBDIRS)

   file( GLOB_RECURSE _tutorials "${_subdir}/*CMakeLists.txt" )

   foreach ( _item  IN LISTS _tutorials)
      get_filename_component(_dir ${_item} DIRECTORY )
      add_subdirectory(${_dir})
   endforeach ()

endforeach ()




#
# Create a target to build all tutorials and an indentical one for
# test_install
#
# !!!!!! ATTENTION: target "tutorial" and "test_install" are the same for
#  the time being until we decide what to do with test_install
#
# add_custom_target (tutorials)
# add_custom_target (test_install)

#
# Here we define a macro to accumulate all the tutorial targets
#
set ( TUTORIAL_TARGETS )

# macro (add_tutorial targets_list)
#    # add_dependencies (tutorials ${targets_list})
#    # add_dependencies (test_install ${targets_list})
# endmacro ()



# # Just write something if
# add_custom_command(TARGET test_install POST_BUILD
#                    COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan
#                    "AMReX installation is working")
