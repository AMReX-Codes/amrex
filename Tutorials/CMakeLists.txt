#
# List of subdirectories to search for CMakeLists.
# For now, we do not include MUI, SDC, SWFFT
#
set( AMREX_TUTORIALS_SUBDIRS Amr Basic ForkJoin)

if (AMReX_ASCENT)
   list(APPEND AMREX_TUTORIALS_SUBDIRS Blueprint)
endif ()

if (AMReX_EB)
   list(APPEND AMREX_TUTORIALS_SUBDIRS EB)
endif ()

if (AMReX_LINEAR_SOLVERS)
   list(APPEND AMREX_TUTORIALS_SUBDIRS LinearSolvers)
endif ()

if (AMReX_PARTICLES)
   list(APPEND AMREX_TUTORIALS_SUBDIRS Particles)
endif ()

if (AMReX_FORTRAN_INTERFACES)
   list(APPEND AMREX_TUTORIALS_SUBDIRS FortranInterface)
endif ()

if (AMReX_CUDA)
   list(APPEND AMREX_TUTORIALS_SUBDIRS GPU)
endif ()

list(TRANSFORM AMREX_TUTORIALS_SUBDIRS PREPEND "${CMAKE_CURRENT_LIST_DIR}/")

#
# Function to setup the tutorials
#
function (setup_tutorial _srcs  _inputs)

   cmake_parse_arguments( "" "HAS_FORTRAN_MODULES"
      "BASE_NAME;RUNTIME_SUBDIR;EXTRA_DEFINITIONS" "" ${ARGN} )

   if (_BASE_NAME)
      set(_base_name ${_BASE_NAME})
   else ()
      string(REGEX REPLACE ".*Tutorials/" "" _base_name ${CMAKE_CURRENT_LIST_DIR})
      string(REPLACE "/" "_" _base_name ${_base_name})
   endif ()

   # Prepend "Tutorial_" to base name
   set(_base_name "Tutorial_${_base_name}")

   if (_RUNTIME_SUBDIR)
      set(_exe_dir ${CMAKE_CURRENT_BINARY_DIR}/${_RUNTIME_SUBDIR})
   else ()
      set(_exe_dir ${CMAKE_CURRENT_BINARY_DIR})
   endif ()

   set( _exe_name  ${_base_name} )

   add_executable( ${_exe_name} )

   target_sources( ${_exe_name} PRIVATE ${${_srcs}} )

   set_target_properties( ${_exe_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${_exe_dir} )

   if (_EXTRA_DEFINITIONS)
      target_compile_definitions(${_exe_name} PRIVATE ${_EXTRA_DEFINITIONS})
   endif ()

   # Find out which include directory is needed
   set(_includes ${${_srcs}})
   list(FILTER _includes INCLUDE REGEX "\\.H$")
   foreach(_item IN LISTS _includes)
      get_filename_component( _include_dir ${_item} DIRECTORY)
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

   if (AMReX_CUDA)
      setup_target_for_cuda_compilation( ${_exe_name} )
   endif ()

   if (${_inputs})
      file( COPY ${${_inputs}} DESTINATION ${_exe_dir} )
   endif ()

endfunction ()

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
