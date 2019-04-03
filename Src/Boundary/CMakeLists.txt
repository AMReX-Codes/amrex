target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(amrex
   PRIVATE 
   AMReX_FabSet.cpp
   AMReX_BndryRegister.cpp
   AMReX_Mask.cpp
   AMReX_MultiMask.cpp
   AMReX_BndryData.cpp
   AMReX_InterpBndryData.cpp
   AMReX_MacBndry.cpp
   AMReX_FabSet.H
   AMReX_BndryRegister.H
   AMReX_Mask.H
   AMReX_MultiMask.H
   AMReX_BndryData.H
   AMReX_BoundCond.H
   AMReX_InterpBndryData.H
   AMReX_LO_BCTYPES.H
   AMReX_MacBndry.H
   AMReX_InterpBndryData_K.H
   AMReX_InterpBndryData_${DIM}D_K.H
   AMReX_LO_UTIL.F90
   AMReX_LOUtil_K.H
   AMReX_YAFluxRegister.H
   AMReX_YAFluxRegister_K.H
   AMReX_YAFluxRegister_${DIM}D_K.H
   AMReX_YAFluxRegister.cpp
   AMReX_lo_bctypes_mod.F90
   )

