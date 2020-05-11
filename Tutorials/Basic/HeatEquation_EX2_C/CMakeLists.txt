#
# This tutorial needs Fortran
#
if ( NOT CMAKE_Fortran_COMPILER_LOADED )
   return()
endif ()

get_filename_component( _base_name ${CMAKE_CURRENT_LIST_DIR} NAME )

set( _test_name ${_base_name} )
set( _exe_name  ${_base_name}.exe )
set( _src_dir   ${CMAKE_CURRENT_LIST_DIR}/Source )
set( _exe_dir   ${CMAKE_CURRENT_LIST_DIR}/Exec )

#
# Create target (executable) for this tutorial
#
add_executable( ${_exe_name} EXCLUDE_FROM_ALL "" )

target_sources( ${_exe_name}
   PRIVATE
   ${_src_dir}/main.cpp
   ${_src_dir}/advance.cpp
   ${_src_dir}/myfunc_F.H
   ${_src_dir}/myfunc.H
   ${_src_dir}/init_phi_${DIM}d.f90
   ${_src_dir}/advance_${DIM}d.f90
   )

target_include_directories( ${_exe_name} PRIVATE  ${_src_dir} )

if (ENABLE_CUDA)
   set_cpp_sources_to_cuda_language( ${_exe_name} )
endif ()

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
file( GLOB_RECURSE _input_files LIST_DIRECTORIES false  ${_exe_dir}/input* )
file( COPY ${_exe_dir}/inputs_${DIM}d DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

#
# Add to the "tutorial" target
#
add_tutorial(${_exe_name})

#
# Clean up
#
unset( _test_name )
unset( _exe_name  )
unset( _src_dir   )
unset( _exe_dir   )
unset( _input_files )
