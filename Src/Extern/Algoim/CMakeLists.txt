target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(amrex
   PRIVATE 
   AMReX_algoim_integrals.H
   AMReX_algoim_integrals.cpp
   )

