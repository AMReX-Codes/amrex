set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")


list ( APPEND ALLHEADERS AMReX_LevelBld.H )
list ( APPEND ALLHEADERS AMReX_Amr.H   AMReX_AmrLevel.H   AMReX_Derive.H   AMReX_StateData.H )
list ( APPEND CXXSRC     AMReX_Amr.cpp AMReX_AmrLevel.cpp AMReX_Derive.cpp AMReX_StateData.cpp )

list ( APPEND ALLHEADERS AMReX_PROB_AMR_F.H )
list ( APPEND ALLHEADERS AMReX_StateDescriptor.H   AMReX_AuxBoundaryData.H   AMReX_Extrapolater.H )
list ( APPEND CXXSRC     AMReX_StateDescriptor.cpp AMReX_AuxBoundaryData.cpp AMReX_Extrapolater.cpp )

list ( APPEND F90SRC     AMReX_extrapolater_${DIM}d.f90)

list ( APPEND F77SRC     AMReX_ARRAYLIM_${DIM}D.F)

#
# This is an option to display arrays in GDB. It is for debugging
# AMReX and currently not supported by cmake.# 
# if(USE_ARRAYVIEW EQUAL TRUE)
#   set(CXX_source_files ${CXX_source_files} AMReX_DatasetClient.cpp)
#   set(CXX_header_files ${CXX_header_files} AMReX_DatasetClient.H)
# endif()

# 
#  Collect sources
# 
set ( ALLSRC  ${CXXSRC} ${F90SRC} ${F77SRC} )

#
# Add absolute path (needed by target_sources)
# 
prepend ( ALLSRC     ${CMAKE_CURRENT_LIST_DIR} )
prepend ( ALLHEADERS ${CMAKE_CURRENT_LIST_DIR} )

target_sources ( amrex PRIVATE ${ALLSRC} ${ALLHEADERS} )

# 
# Accumulate headers 
#
accumulate_headers (ALLHEADERS)   

#
# Add to include dirs
#
accumulate_include_dirs (${CMAKE_CURRENT_LIST_DIR})
