add_amrex_define(AMREX_USE_HYPRE NO_LEGACY NO_1D IF AMReX_HYPRE)

foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
        message(WARNING "HYPRE interfaces are not supported for 1D builds (skipping)")
        continue()
    endif ()

    #
    # This file gets processed if either AMReX_PETSC or AMReX_HYPRE are ON
    # If only AMReX_PETSC is ON, we just need to include this directory
    # and a couple of files to the build.
    #
    target_include_directories(amrex_${D}d
       PUBLIC
       $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

    if (NOT AMReX_HYPRE) # only AMReX_PETSC is active

       target_sources(amrex_${D}d
          PRIVATE
          AMReX_Habec_${D}D_K.H
          AMReX_Habec_K.H
          )

    else ()

       target_sources(amrex_${D}d
          PRIVATE
          AMReX_HypreABecLap.cpp
          AMReX_HypreABecLap.H
          AMReX_HypreABecLap2.cpp
          AMReX_HypreABecLap2.H
          AMReX_HypreABecLap3.cpp
          AMReX_HypreABecLap2.H
          AMReX_Hypre.cpp
          AMReX_Hypre.H
          AMReX_Habec_${D}D_K.H
          AMReX_Habec_K.H
          AMReX_HypreNodeLap.cpp
          AMReX_HypreNodeLap.H
          AMReX_HypreIJIface.cpp
          AMReX_HypreIJIface.H
          AMReX_HypreSolver.H
          AMReX_HypreMLABecLap.H
          AMReX_HypreMLABecLap.cpp
          AMReX_HypreMLABecLap_K.H
          AMReX_HypreMLABecLap_${D}D_K.H
          )

    endif ()
endforeach()
