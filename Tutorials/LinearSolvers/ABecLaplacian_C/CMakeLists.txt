if (DIM EQUAL 1)
   return()
endif ()

#
# This test works for both 2D and 3D builds 
#
add_executable( ABecLaplacian_C.exe  EXCLUDE_FROM_ALL "")

target_sources( ABecLaplacian_C.exe
   PRIVATE
   main.cpp
   MyTest.cpp
   initProb.cpp
   MyTestPlotfile.cpp
   MyTest.H
   initProb_K.H
   )

target_include_directories( ABecLaplacian_C.exe
   PRIVATE
   ${CMAKE_CURRENT_LIST_DIR}
   ) 

set_target_properties( ABecLaplacian_C.exe
   PROPERTIES      
   INCLUDE_DIRECTORIES
   ${CMAKE_CURRENT_BINARY_DIR}/mod_files
   Fortran_MODULE_DIRECTORY
   ${CMAKE_CURRENT_BINARY_DIR}/mod_files
   RUNTIME_OUTPUT_DIRECTORY
   ${CMAKE_CURRENT_BINARY_DIR}
   )

target_link_libraries( ABecLaplacian_C.exe
   PRIVATE
   amrex
   )


# Find input files 
file( GLOB_RECURSE inputs LIST_DIRECTORIES false ${CMAKE_CURRENT_LIST_DIR}/input* )

# Copy input files to corresponding build dir
file( COPY ${inputs} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} ) 

# Add to the "tutorial" target
add_tutorial(ABecLaplacian_C.exe)   

