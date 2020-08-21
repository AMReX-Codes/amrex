#
# This tutorial needs Fortran interfaces
# Disable for cuda
#
if ( NOT ENABLE_FORTRAN_INTERFACES OR ENABLE_CUDA )
   return()
endif ()

get_filename_component( _base_name ${CMAKE_CURRENT_LIST_DIR} NAME )

set( _test_name ${_base_name} )
set( _exe_name  ${_base_name}.exe )

#
# Create target (executable) for this tutorial
#
add_executable( ${_exe_name} EXCLUDE_FROM_ALL "" )

target_sources( ${_exe_name}
   PRIVATE
   fmain.F90
   init_phi.f90
   advance.f90
   )

target_link_libraries( ${_exe_name} amrex )


set_target_properties ( ${_exe_name} PROPERTIES
   INCLUDE_DIRECTORIES
   ${CMAKE_CURRENT_BINARY_DIR}/mod_files
   Fortran_MODULE_DIRECTORY
   ${CMAKE_CURRENT_BINARY_DIR}/mod_files
   )

#
# Find input files and copy them to corresponding build dir
#
file( GLOB_RECURSE _input_files LIST_DIRECTORIES false  input* )
file( COPY ${_input_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

#
# Add to the "tutorial" target
#
add_tutorial(${_exe_name})

#
# Clean up
#
unset( _test_name )
unset( _exe_name  )
unset( _input_files )
