target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)
add_amrex_define(BL_USE_SENSEI_INSITU NO_LEGACY)

target_sources(amrex
   PRIVATE
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
