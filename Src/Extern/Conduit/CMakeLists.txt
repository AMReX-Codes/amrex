target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

# TODO: Particles PR merges another file
target_sources(amrex
   PRIVATE
   AMReX_Conduit_Blueprint.H
   AMReX_Conduit_Blueprint.cpp
   )
