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
   AMReX_FillPatchUtil_I.H
   AMReX_FillPatcher.H
   AMReX_FluxRegister.H
   AMReX_InterpBase.H
   AMReX_InterpBase.cpp
   AMReX_MFInterpolater.H
   AMReX_MFInterpolater.cpp
   AMReX_Interpolater.cpp
   AMReX_TagBox.cpp
   AMReX_AmrMesh.cpp
   AMReX_Interpolater.H
   AMReX_TagBox.H
   AMReX_AmrMesh.H
   AMReX_FluxReg_${AMReX_SPACEDIM}D_C.H
   AMReX_FluxReg_C.H
   AMReX_Interp_C.H
   AMReX_Interp_${AMReX_SPACEDIM}D_C.H
   AMReX_MFInterp_C.H
   AMReX_MFInterp_${AMReX_SPACEDIM}D_C.H
   AMReX_InterpFaceRegister.H
   AMReX_InterpFaceRegister.cpp
   AMReX_InterpFaceReg_C.H
   AMReX_InterpFaceReg_${AMReX_SPACEDIM}D_C.H
   AMReX_AmrCoreFwd.H
   )

if (AMReX_FORTRAN)
   target_sources(amrex
      PRIVATE
      AMReX_FillPatchUtil_F.H
      AMReX_FillPatchUtil_${AMReX_SPACEDIM}d.F90
      )
endif ()

# These files are needed only if Particles are enabled
if ( AMReX_PARTICLES )
  target_sources( amrex PRIVATE AMReX_AmrParGDB.H AMReX_AmrParticles.H )
endif ()
