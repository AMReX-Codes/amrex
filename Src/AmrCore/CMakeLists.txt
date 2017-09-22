set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")


set ( CXXSRC
   AMReX_AmrCore.cpp AMReX_Cluster.cpp AMReX_ErrorList.cpp
   AMReX_FillPatchUtil.cpp AMReX_FluxRegister.cpp
   AMReX_Interpolater.cpp AMReX_TagBox.cpp AMReX_AmrMesh.cpp AMReX_YAFluxRegister.cpp )

set ( F77SRC AMReX_FLUXREG_${BL_SPACEDIM}D.F AMReX_INTERP_${BL_SPACEDIM}D.F )

set ( F90SRC AMReX_FillPatchUtil_${BL_SPACEDIM}d.F90 AMReX_YAFluxRegister_$(BL_SPACEDIM)d.F90 AMReX_YAFluxRegister_nd.F90 )

set ( ALLHEADERS
   AMReX_AmrCore.H AMReX_Cluster.H AMReX_ErrorList.H
   AMReX_FillPatchUtil.H AMReX_FluxRegister.H AMReX_Interpolater.H
   AMReX_TagBox.H AMReX_AmrMesh.H AMReX_INTERP_F.H
   AMReX_FillPatchUtil_F.H AMReX_FLUXREG_F.H
   AMReX_FillPatchUtil_F.H AMReX_YAFluxRegister_F.H AMReX_YAFluxRegister.H)


# These files are neeeded only if Particles are enabled
if (ENABLE_PARTICLES EQUAL 1)
  list ( APPEND ALLHEADERS AMReX_AmrParGDB.H AMReX_AmrParticles.H )
endif ()


# Collect sources
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
