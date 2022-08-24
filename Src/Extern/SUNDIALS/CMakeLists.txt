target_include_directories( amrex
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( amrex
   PRIVATE
   AMReX_Sundials.H
   AMReX_Sundials.cpp
   AMReX_NVector_MultiFab.cpp
   AMReX_NVector_MultiFab.H
   AMReX_SUNMemory.cpp
   AMReX_SUNMemory.H
   AMReX_SundialsIntegrator.H
   )
