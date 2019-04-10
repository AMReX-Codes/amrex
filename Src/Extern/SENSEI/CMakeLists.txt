find_package(SENSEI REQUIRED)

target_link_libraries( amrex PUBLIC sensei)   
target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)
target_compile_definitions(amrex PUBLIC $<BUILD_INTERFACE:BL_USE_SENSEI_INSITU>)

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
