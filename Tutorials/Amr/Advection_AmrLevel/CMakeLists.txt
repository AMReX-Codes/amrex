#
# This test works for both 2D and 3D builds
# There are two different tutorials here:
#
# 1) Single Vortex
# 2) Uniform Velocity
#

set (SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/Source )

#
# Sources common to both cases
# 
scan_for_sources (F90SRC F77SRC CXXSRC ALLHEADERS
   ROOT ${SRC_DIR}/Src_${DIM}d )
         
scan_for_sources (TMP1 TMP2 TMP3 TMP4
   ROOT ${SRC_DIR}/Src_nd )

list (APPEND F90SRC     ${TMP1} )
list (APPEND F77SRC     ${TMP2} )
list (APPEND CXXSRC     ${TMP3} ${SRC_DIR}/AmrLevelAdv.cpp  ${SRC_DIR}/LevelBldAdv.cpp
   ${SRC_DIR}/main.cpp )
list (APPEND ALLHEADERS ${TMP4} ${SRC_DIR}/Adv_F.H ${SRC_DIR}/AmrLevelAdv.H)


#
# Single Vortex tutorial ------------------------------------------------------
# 
set ( EXENAME  "Advection_AmrLevel_SV.exe")

set (EXE_DIR ${CMAKE_CURRENT_LIST_DIR}/Exec/SingleVortex)

list ( APPEND F90SRC_SV ${F90SRC} ${EXE_DIR}/face_velocity_${DIM}d.f90
   ${EXE_DIR}/Prob.f90 )

find_include_paths (TUTORIAL_INCLUDE_PATHS ROOT ${CMAKE_CURRENT_LIST_DIR})

#
# Find input files 
#
file ( GLOB_RECURSE inputs LIST_DIRECTORIES false 
   ${EXE_DIR}/input* ${EXE_DIR}/probin)

#
# Copy input files to corresponding build dir
#
file ( COPY ${inputs} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/SingleVortex ) 

#
# Create executable 
# 
add_executable ( ${EXENAME} EXCLUDE_FROM_ALL
   ${F90SRC_SV} ${F77SRC} ${CXXSRC} ${ALLHEADERS} )

target_link_libraries ( ${EXENAME} amrex ${AMREX_EXTRA_CXX_LINK_LINE} )

# Fortran modules for the target will not be added to the amrex mod_files folder 
set ( MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files )

# Always put MOD_DIR before AMReX includes paths to avoid conflicts
# between amrex fortran modules names and tutorial fortran modules names.  
set_target_properties ( ${EXENAME} PROPERTIES
   INCLUDE_DIRECTORIES "${MOD_DIR};${AMREX_INCLUDE_PATHS};${TUTORIAL_INCLUDE_PATHS}"
   Fortran_MODULE_DIRECTORY ${MOD_DIR}
   RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/SingleVortex )

add_tutorial (${EXENAME})   

#
# Uniform Velocity tutorial ------------------------------------------------------
# 
set ( EXENAME  "Advection_AmrLevel_UV.exe")

set (EXE_DIR ${CMAKE_CURRENT_LIST_DIR}/Exec/UniformVelocity)

list ( APPEND F90SRC_UV ${F90SRC} ${EXE_DIR}/face_velocity_${DIM}d.f90
   ${EXE_DIR}/Prob.f90 ${EXE_DIR}/probdata.f90)

find_include_paths (TUTORIAL_INCLUDE_PATHS ROOT ${CMAKE_CURRENT_LIST_DIR})

#
# Find input files 
#
file ( GLOB_RECURSE inputs LIST_DIRECTORIES false 
   ${EXE_DIR}/input* ${EXE_DIR}/probin)

#
# Copy input files to corresponding build dir
#
file ( COPY ${inputs} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/UniformVelocity ) 

#
# Create executable 
# 
add_executable ( ${EXENAME} EXCLUDE_FROM_ALL
   ${F90SRC_UV} ${F77SRC} ${CXXSRC} ${ALLHEADERS} )

target_link_libraries ( ${EXENAME} amrex ${AMREX_EXTRA_CXX_LINK_LINE} )


# Fortran modules for the target will not be added to the amrex mod_files folder 
set ( MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files )

# Always put MOD_DIR before AMReX includes paths to avoid conflicts
# between amrex fortran modules names and tutorial fortran modules names.  
set_target_properties ( ${EXENAME} PROPERTIES
   INCLUDE_DIRECTORIES "${MOD_DIR};${AMREX_INCLUDE_PATHS};${TUTORIAL_INCLUDE_PATHS}"
   Fortran_MODULE_DIRECTORY ${MOD_DIR}
   RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/UniformVelocity )

add_tutorial (${EXENAME})
