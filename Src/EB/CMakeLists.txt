set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")

list ( APPEND ALLHEADERS AMReX_EBAmrUtil_F.H     AMReX_EBCellFlag.H        AMReX_EBFArrayBox.H )
list ( APPEND ALLHEADERS AMReX_EBInterp_F.H      AMReX_EBMultiFabUtil.H    AMReX_MultiCutFab.H )
list ( APPEND ALLHEADERS AMReX_EBAmrUtil.H       AMReX_EBDataCollection.H  AMReX_EBFluxRegister_F.H )
list ( APPEND ALLHEADERS AMReX_EBInterpolater.H  AMReX_EBSupport.H         AMReX_EBCellFlag_F.H )
list ( APPEND ALLHEADERS AMReX_EBFabFactory.H    AMReX_EBFluxRegister.H    AMReX_EBMultiFabUtil_F.H  AMReX_EBTower.H )
list ( APPEND ALLHEADERS AMReX_EB_F.H            AMReX_EB_levelset.H )

list ( APPEND CXXSRC AMReX_EBAmrUtil.cpp       AMReX_EBDataCollection.cpp  AMReX_EBFArrayBox.cpp )
list ( APPEND CXXSRC AMReX_EBInterpolater.cpp  AMReX_EBTower.cpp )
list ( APPEND CXXSRC AMReX_EBCellFlag.cpp      AMReX_EBFabFactory.cpp      AMReX_EBFluxRegister.cpp   )
list ( APPEND CXXSRC AMReX_EBMultiFabUtil.cpp  AMReX_MultiCutFab.cpp )
list ( APPEND CXXSRC AMReX_EB_levelset.cpp )

list ( APPEND F90SRC AMReX_EBFluxRegister_${DIM}d.F90 )
list ( APPEND F90SRC AMReX_EBFluxRegister_nd.F90 )
list ( APPEND F90SRC AMReX_EBMultiFabUtil_${DIM}d.F90 )
list ( APPEND F90SRC AMReX_EBMultiFabUtil_nd.F90 )
list ( APPEND F90SRC AMReX_EBAmrUtil_nd.F90      )
list ( APPEND F90SRC AMReX_ebinterp_${DIM}d.F90   )
list ( APPEND F90SRC AMReX_ebcellflag_mod.F90    )
list ( APPEND F90SRC AMReX_compute_normals.F90 AMReX_EB_geometry.F90 AMReX_EB_levelset_F.F90 )

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
