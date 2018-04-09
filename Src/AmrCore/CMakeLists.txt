set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")



list ( APPEND CXXSRC      AMReX_AmrCore.cpp AMReX_Cluster.cpp AMReX_ErrorList.cpp ) 
list ( APPEND ALLHEADERS  AMReX_AmrCore.H   AMReX_Cluster.H   AMReX_ErrorList.H )


list ( APPEND CXXSRC      AMReX_FillPatchUtil.cpp AMReX_FluxRegister.cpp )
list ( APPEND ALLHEADERS  AMReX_FillPatchUtil.H   AMReX_FluxRegister.H )


list ( APPEND CXXSRC      AMReX_Interpolater.cpp AMReX_TagBox.cpp AMReX_AmrMesh.cpp )
list ( APPEND ALLHEADERS  AMReX_Interpolater.H   AMReX_TagBox.H   AMReX_AmrMesh.H   )

list ( APPEND F90SRC      AMReX_FLUXREG_${DIM}D.F90  AMReX_FLUXREG_nd.F90  AMReX_INTERP_${DIM}D.F90 )
list ( APPEND ALLHEADERS  AMReX_FLUXREG_F.H        AMReX_INTERP_F.H )

list ( APPEND F90SRC      AMReX_FillPatchUtil_${DIM}d.F90 )
list ( APPEND ALLHEADERS  AMReX_FillPatchUtil_F.H )


# These files are neeeded only if Particles are enabled
if ( ENABLE_PARTICLES )
  list ( APPEND ALLHEADERS AMReX_AmrParGDB.H AMReX_AmrParticles.H )
endif ()

# 
# Collect sources
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
