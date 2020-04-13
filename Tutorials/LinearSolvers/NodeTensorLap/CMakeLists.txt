if (DIM EQUAL 1)
   return()
endif ()

#
# Use directory name as exe name
#
get_filename_component(exe_name ${CMAKE_CURRENT_LIST_DIR} NAME)
set(exe_name "${exe_name}.exe")

#
# This test works for both 2D and 3D builds
#
add_executable( ${exe_name} EXCLUDE_FROM_ALL "")

target_sources( ${exe_name}
   PRIVATE
   main.cpp
   MyTest.cpp
   MyTest.H
   MyTestPlotfile.cpp
   )

if (ENABLE_CUDA)
   set_cpp_sources_to_cuda_language( ${exe_name} )
endif ()

target_include_directories(${exe_name} PRIVATE ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries( ${exe_name} PRIVATE amrex )

# Find input files
file( GLOB_RECURSE inputs LIST_DIRECTORIES false ${CMAKE_CURRENT_LIST_DIR}/input* )

# Copy input files to corresponding build dir
file( COPY ${inputs} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

# Add to the "tutorial" target
add_tutorial(${exe_name})
