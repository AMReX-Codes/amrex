if (DIM EQUAL 1)
   message(FATAL_ERROR "HYPRE interfaces are not supported for 1D builds")
endif ()

find_package(HYPRE 2.15 REQUIRED)

target_compile_definitions( amrex
   PUBLIC
   $<BUILD_INTERFACE:AMREX_USE_HYPRE>)

target_include_directories( amrex
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

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
   AMReX_HypreABec_F.H
   AMReX_HypreNodeLap.cpp
   AMReX_HypreNodeLap.H
   AMReX_HABEC_${DIM}D.F90
   AMReX_Hypre_fort_mod.F90
   )

target_link_libraries( amrex PUBLIC HYPRE )
