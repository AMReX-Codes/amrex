target_include_directories( amrex
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

add_amrex_define(AMREX_USE_BITTREE NO_LEGACY)

target_sources( amrex
   PRIVATE
   AMReX_Bittree.H
   AMReX_Bittree.cpp
   )
