foreach(D IN LISTS AMReX_SPACEDIM)
    #
    # Sources in subdirectory MLMG
    #
    target_include_directories(amrex_${D}d PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)
    target_include_directories(amrex_${D}d PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/MLMG>)

    target_sources(amrex_${D}d
       PRIVATE
       MLMG/AMReX_MLMG.H
       MLMG/AMReX_MLMG.cpp
       MLMG/AMReX_MLMG_K.H
       MLMG/AMReX_MLMG_${D}D_K.H
       MLMG/AMReX_MLMGBndry.H
       MLMG/AMReX_MLLinOp.H
       MLMG/AMReX_MLLinOp_K.H
       MLMG/AMReX_MLCellLinOp.H
       MLMG/AMReX_MLNodeLinOp.H
       MLMG/AMReX_MLNodeLinOp_K.H
       MLMG/AMReX_MLNodeLinOp_${D}D_K.H
       MLMG/AMReX_MLNodeLinOp.cpp
       MLMG/AMReX_MLCellABecLap.H
       MLMG/AMReX_MLCellABecLap_K.H
       MLMG/AMReX_MLCellABecLap_${D}D_K.H
       MLMG/AMReX_MLCGSolver.H
       MLMG/AMReX_PCGSolver.H
       MLMG/AMReX_MLABecLaplacian.H
       MLMG/AMReX_MLABecLap_K.H
       MLMG/AMReX_MLABecLap_${D}D_K.H
       MLMG/AMReX_MLALaplacian.H
       MLMG/AMReX_MLALap_K.H
       MLMG/AMReX_MLALap_${D}D_K.H
       MLMG/AMReX_MLPoisson.H
       MLMG/AMReX_MLPoisson_K.H
       MLMG/AMReX_MLPoisson_${D}D_K.H
       AMReX_GMRES.H
       AMReX_GMRES_MLMG.H
       AMReX_GMRES_MV.H
       AMReX_Smoother_MV.H
       AMReX_Algebra.H
       AMReX_AlgPartition.H
       AMReX_AlgPartition.cpp
       AMReX_AlgVector.H
       AMReX_SpMatrix.H
       AMReX_SpMV.H
       )

    if (D EQUAL 3)
       target_sources(amrex_${D}d
          PRIVATE
          MLMG/AMReX_MLMG_2D_K.H
          MLMG/AMReX_MLPoisson_2D_K.H
          MLMG/AMReX_MLALap_2D_K.H
          )
    endif ()

    if (AMReX_LINEAR_SOLVERS_EM)
       target_sources(amrex_${D}d
          PRIVATE
          MLMG/AMReX_MLCurlCurl.H
          MLMG/AMReX_MLCurlCurl.cpp
          MLMG/AMReX_MLCurlCurl_K.H
          )

       target_sources(amrex_${D}d
          PRIVATE
          MLMG/AMReX_MLEBNodeFDLaplacian.H
          MLMG/AMReX_MLEBNodeFDLaplacian.cpp
          MLMG/AMReX_MLEBNodeFDLap_K.H
          MLMG/AMReX_MLEBNodeFDLap_${D}D_K.H
          MLMG/AMReX_MLNodeTensorLaplacian.H
          MLMG/AMReX_MLNodeTensorLaplacian.cpp
          MLMG/AMReX_MLNodeTensorLap_K.H
          MLMG/AMReX_MLNodeTensorLap_${D}D_K.H
       )
    endif ()

    if (AMReX_LINEAR_SOLVERS_INCFLO)
       target_sources(amrex_${D}d
          PRIVATE
          MLMG/AMReX_MLNodeABecLaplacian.H
          MLMG/AMReX_MLNodeABecLaplacian.cpp
          MLMG/AMReX_MLNodeABecLap_K.H
          MLMG/AMReX_MLNodeABecLap_${D}D_K.H
          MLMG/AMReX_MLNodeLaplacian.H
          MLMG/AMReX_MLNodeLaplacian.cpp
          MLMG/AMReX_MLNodeLaplacian_sync.cpp
          MLMG/AMReX_MLNodeLaplacian_sten.cpp
          MLMG/AMReX_MLNodeLaplacian_misc.cpp
          MLMG/AMReX_MLNodeLap_K.H
          MLMG/AMReX_MLNodeLap_${D}D_K.H
          MLMG/AMReX_MLTensorOp.H
          MLMG/AMReX_MLTensorOp.cpp
          MLMG/AMReX_MLTensorOp_grad.cpp
          MLMG/AMReX_MLTensor_K.H
          MLMG/AMReX_MLTensor_${D}D_K.H
       )
    endif ()

    if (AMReX_EB AND NOT D EQUAL 1)
       target_sources(amrex_${D}d
          PRIVATE
          MLMG/AMReX_MLEBABecLap.H
          MLMG/AMReX_MLEBABecLap.cpp
          MLMG/AMReX_MLEBABecLap_F.cpp
          MLMG/AMReX_MLEBABecLap_K.H
          MLMG/AMReX_MLEBABecLap_${D}D_K.H
          )

       if (AMReX_LINEAR_SOLVERS_INCFLO)
          target_sources(amrex_${D}d
             PRIVATE
             MLMG/AMReX_MLNodeLaplacian_eb.cpp
             MLMG/AMReX_MLEBTensorOp.H
             MLMG/AMReX_MLEBTensorOp.cpp
             MLMG/AMReX_MLEBTensorOp_bc.cpp
             MLMG/AMReX_MLEBTensor_K.H
             MLMG/AMReX_MLEBTensor_${D}D_K.H
             )
       endif ()
    endif ()

    if (AMReX_FORTRAN)
       target_sources(amrex_${D}d
          PRIVATE
          MLMG/AMReX_MLLinOp_F.H
          MLMG/AMReX_MLLinOp_nd.F90
          )
    endif ()

    if (AMReX_HYPRE AND AMReX_LINEAR_SOLVERS_INCFLO)
       target_sources(amrex_${D}d
          PRIVATE
          MLMG/AMReX_MLNodeLaplacian_hypre.cpp
          )
    endif ()

    if (D EQUAL 3)

       target_include_directories(amrex_${D}d PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/OpenBC>)

       target_sources(amrex_${D}d
          PRIVATE
          OpenBC/AMReX_OpenBC.H
          OpenBC/AMReX_OpenBC_K.H
          OpenBC/AMReX_OpenBC.cpp
          )
    endif ()
endforeach()
