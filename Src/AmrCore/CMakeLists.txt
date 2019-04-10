target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(amrex
   PRIVATE 
   AMReX_AmrCore.cpp
   AMReX_Cluster.cpp
   AMReX_ErrorList.cpp 
   AMReX_AmrCore.H
   AMReX_Cluster.H
   AMReX_ErrorList.H
   AMReX_FillPatchUtil.cpp
   AMReX_FluxRegister.cpp
   AMReX_FillPatchUtil.H
   AMReX_FluxRegister.H
   AMReX_Interpolater.cpp
   AMReX_TagBox.cpp
   AMReX_AmrMesh.cpp
   AMReX_Interpolater.H
   AMReX_TagBox.H
   AMReX_AmrMesh.H 
   AMReX_FluxReg_${DIM}D_C.H
   AMReX_FluxReg_C.H
   AMReX_FLUXREG_nd.F90
   AMReX_INTERP_${DIM}D.F90
   AMReX_FLUXREG_F.H
   AMReX_INTERP_F.H
   AMReX_Interp_C.H
   AMReX_Interp_${DIM}D_C.H
   AMReX_FillPatchUtil_${DIM}d.F90
   AMReX_FillPatchUtil_F.H )


# These files are neeeded only if Particles are enabled
if ( ENABLE_PARTICLES )
  target_sources( amrex PRIVATE AMReX_AmrParGDB.H AMReX_AmrParticles.H )
endif ()
