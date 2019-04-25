target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(amrex
   PRIVATE 
   AMReX_LevelBld.H
   AMReX_Amr.H
   AMReX_AmrLevel.H
   AMReX_Derive.H
   AMReX_StateData.H
   AMReX_Amr.cpp
   AMReX_AmrLevel.cpp
   AMReX_Derive.cpp
   AMReX_StateData.cpp
   AMReX_PROB_AMR_F.H
   AMReX_StateDescriptor.H
   AMReX_AuxBoundaryData.H
   AMReX_Extrapolater.H
   AMReX_StateDescriptor.cpp
   AMReX_AuxBoundaryData.cpp
   AMReX_Extrapolater.cpp
   AMReX_extrapolater_${DIM}d.f90 )
