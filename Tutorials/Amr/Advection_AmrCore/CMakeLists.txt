if (DIM EQUAL 1)
   return()
endif ()

#
# This test works for both 2D and 3D builds 
#
set ( EXENAME  "Advection_AmrCore.exe" )
set ( SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/Source )
set ( EXE_DIR ${CMAKE_CURRENT_LIST_DIR}/Exec/SingleVortex )

#
# Create target (executable) for this tutorial
# 
add_executable ( ${EXENAME} EXCLUDE_FROM_ALL "" )

#
# Set target (executable) properties 
#
find_all_sources ( SRC1 INC1 ROOT ${SRC_DIR} )
find_all_sources ( SRC2 INC2 ROOT ${SRC_DIR}/Src_${DIM}d )
find_all_sources ( SRC3 INC3 ROOT ${SRC_DIR}/Src_nd )

target_sources ( ${EXENAME}
   PRIVATE
   ${EXE_DIR}/Prob.f90 ${EXE_DIR}/face_velocity_${DIM}d.f90
   ${SRC1} ${SRC2} ${SRC3} )

set_target_properties ( ${EXENAME} PROPERTIES      
   INCLUDE_DIRECTORIES
   "${INC1};${INC2};${INC3};${CMAKE_CURRENT_BINARY_DIR}/mod_files"
   Fortran_MODULE_DIRECTORY
   ${CMAKE_CURRENT_BINARY_DIR}/mod_files
   RUNTIME_OUTPUT_DIRECTORY
   ${CMAKE_CURRENT_BINARY_DIR} )

target_link_libraries ( ${EXENAME} amrex ${AMREX_EXTRA_Fortran_LINK_LINE} )

#
# Find input files 
#
file ( GLOB_RECURSE inputs LIST_DIRECTORIES false  ${EXE_DIR}/input* )

#
# Copy input files to corresponding build dir
#
file ( COPY ${inputs} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} ) 

#
# Add to the "tutorial" target
# 
add_tutorial (${EXENAME})   

