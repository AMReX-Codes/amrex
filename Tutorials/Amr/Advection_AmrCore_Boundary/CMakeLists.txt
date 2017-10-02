#
# This test works for both 2D and 3D builds
#
set ( EXENAME  "Advection_AmrCore_Boundary.exe")
set (SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/Source )
set (EXE_DIR ${CMAKE_CURRENT_LIST_DIR}/Exec/SingleVortex)

scan_for_sources (F90SRC F77SRC CXXSRC ALLHEADERS
   ROOT ${SRC_DIR}/Src_${BL_SPACEDIM}d )
         
scan_for_sources (TMP1 TMP2 TMP3 TMP4
   ROOT ${SRC_DIR}/Src_nd )

list (APPEND F90SRC     ${TMP1} ${EXE_DIR}/Prob.f90 ${EXE_DIR}/face_velocity_${BL_SPACEDIM}d.f90)
list (APPEND F77SRC     ${TMP2})
list (APPEND CXXSRC     ${TMP3} ${SRC_DIR}/AmrCoreAdv.cpp  ${SRC_DIR}/main.cpp)
list (APPEND ALLHEADERS ${TMP4} ${SRC_DIR}/AmrCoreAdv_F.H ${SRC_DIR}/AmrCoreAdv.H
   ${SRC_DIR}/AmrCoreAdvPhysBC.H)

find_include_paths (TUTORIAL_INCLUDE_PATHS ROOT ${CMAKE_CURRENT_LIST_DIR})

#
# Find input files 
#
file ( GLOB_RECURSE inputs LIST_DIRECTORIES false 
   ${EXE_DIR}/input* )

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
   INCLUDE_DIRECTORIES "${AMREX_INCLUDE_PATHS};${TUTORIAL_INCLUDE_PATHS}"
   RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )

add_tutorial (${EXENAME})   

