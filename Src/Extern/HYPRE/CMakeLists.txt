add_amrex_define(AMREX_USE_HYPRE NO_LEGACY IF AMReX_HYPRE)

foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
        continue()
    endif ()

    #
    # This file gets processed if either AMReX_PETSC or AMReX_HYPRE are ON.
    # The matrix kernels shared with PETSc live in LinearSolvers/MLMG.
    #
    target_include_directories(amrex_${D}d
       PUBLIC
       $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

    if (AMReX_HYPRE)

       target_sources(amrex_${D}d
          PRIVATE
          AMReX_HypreABecLap.cpp
          AMReX_HypreABecLap.H
          AMReX_HypreABecLap2.cpp
          AMReX_HypreABecLap2.H
          AMReX_HypreABecLap3.cpp
          AMReX_HypreABecLap3.H
          AMReX_Hypre.cpp
          AMReX_Hypre.H
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
