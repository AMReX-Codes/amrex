#
# This test requires a 2D build
# 
if ( NOT (${BL_SPACEDIM} EQUAL 2) )
   return ()
endif ()

set ( EXENAME  "Advection_F.exe")


#
#  We use this very bad method to find sources
#  automagically
#
scan_for_sources (F90SRC F77SRC CXXSRC ALLHEADERS)

#
# Find input files 
#
file ( GLOB_RECURSE inputs LIST_DIRECTORIES false 
   ${CMAKE_CURRENT_LIST_DIR}/input* )

#
# Copy input files to corresponding build dir
#
file ( COPY ${inputs} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} ) 

#
# Create executable 
# 
add_executable ( ${EXENAME} EXCLUDE_FROM_ALL
   ${F90SRC} ${F77SRC} ${CXXSRC} ${ALLHEADERS} )

target_link_libraries ( ${EXENAME} amrex ${AMREX_EXTRA_Fortran_LINK_LINE} )

set_target_properties ( ${EXENAME} PROPERTIES
   INCLUDE_DIRECTORIES "${AMREX_INCLUDE_PATHS}"
   RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )

add_tutorial (${EXENAME})   

