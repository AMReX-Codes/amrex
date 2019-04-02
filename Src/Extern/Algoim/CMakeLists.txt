target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_compile_definitions(amrex PUBLIC $<BUILD_INTERFACE:USE_ALGOIM>)

target_sources(amrex
   PRIVATE 
   AMReX_algoim_integrals.H
   AMReX_algoim_integrals.cpp
   )

