target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)
add_amrex_define(BL_USE_SENSEI_INSITU NO_LEGACY)

# add option for automatic instrumentation of AMR
add_amrex_define(AMREX_USE_SENSEI_AUTO IF AMReX_SENSEI_AUTO NO_LEGACY)

set (amrex_sensei_sources
   AMReX_AmrMeshDataAdaptor.H
   AMReX_AmrDataAdaptor.H
   AMReX_InSituBridge.H
   AMReX_AmrMeshInSituBridge.H
   AMReX_AmrInSituBridge.H
   AMReX_InSituUtils.H
   AMReX_AmrMeshDataAdaptor.cpp
   AMReX_AmrDataAdaptor.cpp
   AMReX_InSituBridge.cpp
   AMReX_AmrMeshInSituBridge.cpp
   AMReX_AmrInSituBridge.cpp
   AMReX_InSituUtils.cpp
)

if (AMReX_PARTICLES)
   list (APPEND amrex_sensei_sources
      AMReX_AmrParticleDataAdaptor.H
      AMReX_AmrParticleDataAdaptorI.H
      AMReX_AmrParticleMeshDataAdaptor.H
      AMReX_AmrParticleMeshDataAdaptorI.H
      AMReX_AmrParticleInSituBridge.H
      AMReX_AmrParticleInSituBridge.cpp
      AMReX_AmrParticleMeshInSituBridge.H
      AMReX_AmrParticleMeshInSituBridge.cpp
   )
endif ()

target_sources(amrex
   PRIVATE
   ${amrex_sensei_sources}
)
