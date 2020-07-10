#
# This test requires a 2D build
# 
if ( NOT (${DIM} EQUAL 2) OR (NOT ENABLE_FORTRAN_INTERFACES) )
   return ()
endif ()

set ( EXENAME  "Advection_octreeF_SingleVortex.exe" )
set ( EXE_DIR ${CMAKE_CURRENT_LIST_DIR}/Exec/SingleVortex)
set ( SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/Source )


add_executable ( ${EXENAME} EXCLUDE_FROM_ALL "" )

#
# Find sources and include paths
#
find_all_sources ( SRC1 INC1 ROOT ${SRC_DIR} )
find_all_sources ( SRC2 INC2 ROOT ${SRC_DIR}/Src_${DIM}d )

#
# Add sources to target
# 
target_sources ( ${EXENAME}
   PRIVATE
   ${EXE_DIR}/Prob.f90 ${EXE_DIR}/face_velocity_${DIM}d.F90
   ${SRC1} ${SRC2} ${SRC3} )

if ("${INC1} ${INC2} ${INC3}")
   target_include_directories ( ${EXENAME} PRIVATE ${INC1} ${INC2} ${INC3} )
endif ()

set_target_properties ( ${EXENAME} PROPERTIES      
   Fortran_MODULE_DIRECTORY
   ${CMAKE_CURRENT_BINARY_DIR}/mod_files
   RUNTIME_OUTPUT_DIRECTORY
   ${CMAKE_CURRENT_BINARY_DIR}/SingleVortex
   LINKER_LANGUAGE
   Fortran
   )

target_link_libraries (${EXENAME} amrex)

#
# Find input files 
#
file ( GLOB_RECURSE inputs LIST_DIRECTORIES false 
   ${CMAKE_CURRENT_LIST_DIR}/input* )

#
# Copy input files to corresponding build dir
#
file ( COPY ${inputs} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/SingleVortex ) 


add_tutorial (${EXENAME})   

