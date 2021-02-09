#
# This file gets processed if either AMReX_PETSC or AMReX_HYPRE are ON
# If only AMReX_PETSC is ON, we just need to include this directory
# and a couple of files to the build.
#
target_include_directories( amrex
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

if (NOT AMReX_HYPRE) # only AMReX_PETSC is active

   target_sources( amrex
      PRIVATE
      AMReX_Habec_${AMReX_SPACEDIM}D_K.H
      AMReX_Habec_K.H
      )

else ()

   if (AMReX_SPACEDIM EQUAL 1)
      message(FATAL_ERROR "HYPRE interfaces are not supported for 1D builds")
   endif ()

   add_amrex_define(AMREX_USE_HYPRE NO_LEGACY)

   target_sources( amrex
      PRIVATE
      AMReX_HypreABecLap.cpp
      AMReX_HypreABecLap.H
      AMReX_HypreABecLap2.cpp
      AMReX_HypreABecLap2.H
      AMReX_HypreABecLap3.cpp
      AMReX_HypreABecLap2.H
      AMReX_Hypre.cpp
      AMReX_Hypre.H
      AMReX_Habec_${AMReX_SPACEDIM}D_K.H
      AMReX_Habec_K.H
      AMReX_HypreNodeLap.cpp
      AMReX_HypreNodeLap.H
      AMReX_HypreIJIface.cpp
      AMReX_HypreIJIface.H
      )

endif ()
