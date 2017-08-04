set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")

set (CXXSRC
   AMReX_Amr.cpp AMReX_AmrLevel.cpp AMReX_AuxBoundaryData.cpp AMReX_Derive.cpp
   AMReX_Extrapolater.cpp AMReX_StateData.cpp AMReX_StateDescriptor.cpp )

set (F90SRC AMReX_extrapolater_${BL_SPACEDIM}d.f90)

set (F77SRC AMReX_ARRAYLIM_${BL_SPACEDIM}D.F)

set (ALLHEADERS
   AMReX_Amr.H AMReX_AmrLevel.H AMReX_AuxBoundaryData.H AMReX_Derive.H
   AMReX_Extrapolater.H AMReX_LevelBld.H AMReX_PROB_AMR_F.H AMReX_StateData.H
   AMReX_StateDescriptor.H)

#
# This is an option to display arrays in GDB. It is for debugging
# AMReX and currently not supported by cmake.# 
# if(USE_ARRAYVIEW EQUAL TRUE)
#   set(CXX_source_files ${CXX_source_files} AMReX_DatasetClient.cpp)
#   set(CXX_header_files ${CXX_header_files} AMReX_DatasetClient.H)
# endif()


#  Collect sources
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
