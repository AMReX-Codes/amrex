add_sources ( AMReX_AmrCore.cpp AMReX_Cluster.cpp AMReX_ErrorList.cpp ) 
add_sources ( AMReX_AmrCore.H   AMReX_Cluster.H   AMReX_ErrorList.H )

add_sources ( AMReX_FillPatchUtil.cpp AMReX_FluxRegister.cpp )
add_sources ( AMReX_FillPatchUtil.H   AMReX_FluxRegister.H )

add_sources ( AMReX_Interpolater.cpp AMReX_TagBox.cpp AMReX_AmrMesh.cpp )
add_sources ( AMReX_Interpolater.H   AMReX_TagBox.H   AMReX_AmrMesh.H  )

add_sources ( AMReX_FLUXREG_${DIM}D.F90  AMReX_FLUXREG_nd.F90  AMReX_INTERP_${DIM}D.F90 )
add_sources ( AMReX_FLUXREG_F.H        AMReX_INTERP_F.H )

add_sources ( AMReX_Interp_C.H AMReX_Interp_${DIM}D_C.H )

add_sources ( AMReX_FillPatchUtil_${DIM}d.F90 )
add_sources ( AMReX_FillPatchUtil_F.H )


# These files are neeeded only if Particles are enabled
if ( ENABLE_PARTICLES )
  add_sources ( AMReX_AmrParGDB.H AMReX_AmrParticles.H )
endif ()
