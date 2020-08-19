#
# This file gets processed if either ENABLE_PETSC or ENABLE_HYPRE are ON
# If only ENABLE_PETSC is ON, we just need to include this directory
# and a couple of files to the build.
#
target_include_directories( amrex
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

if (NOT ENABLE_HYPRE) # only ENABLE_PETSC is active

   target_sources( amrex
      PRIVATE
      AMReX_Habec_${DIM}D_K.H
      AMReX_Habec_K.H
      )

else ()

   if (DIM EQUAL 1)
      message(FATAL_ERROR "HYPRE interfaces are not supported for 1D builds")
   endif ()

   target_compile_definitions( amrex
      PUBLIC
      $<BUILD_INTERFACE:AMREX_USE_HYPRE>)

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
      AMReX_Habec_${DIM}D_K.H
      AMReX_Habec_K.H
      AMReX_HypreNodeLap.cpp
      AMReX_HypreNodeLap.H
      )

endif ()
