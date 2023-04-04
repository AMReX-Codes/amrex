#
# Sources in subdirectory MLMG
#
target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/MLMG>)

target_sources(amrex
   PRIVATE
   MLMG/AMReX_MLMG.H
   MLMG/AMReX_MLMG.cpp
   MLMG/AMReX_MLMG_K.H
   MLMG/AMReX_MLMG_${AMReX_SPACEDIM}D_K.H
   MLMG/AMReX_MLMGBndry.H
   MLMG/AMReX_MLLinOp.H
   MLMG/AMReX_MLLinOp_K.H
   MLMG/AMReX_MLCellLinOp.H
   MLMG/AMReX_MLNodeLinOp.H
   MLMG/AMReX_MLNodeLinOp.cpp
   MLMG/AMReX_MLCellABecLap.H
   MLMG/AMReX_MLCellABecLap_K.H
   MLMG/AMReX_MLCellABecLap_${AMReX_SPACEDIM}D_K.H
   MLMG/AMReX_MLCGSolver.H
   MLMG/AMReX_MLABecLaplacian.H
   MLMG/AMReX_MLABecLap_K.H
   MLMG/AMReX_MLABecLap_${AMReX_SPACEDIM}D_K.H
   MLMG/AMReX_MLALaplacian.H
   MLMG/AMReX_MLALap_K.H
   MLMG/AMReX_MLALap_${AMReX_SPACEDIM}D_K.H
   MLMG/AMReX_MLPoisson.H
   MLMG/AMReX_MLPoisson_K.H
   MLMG/AMReX_MLPoisson_${AMReX_SPACEDIM}D_K.H
   MLMG/AMReX_MLNodeLaplacian.H
   MLMG/AMReX_MLNodeLaplacian.cpp
   MLMG/AMReX_MLNodeLaplacian_sync.cpp
   MLMG/AMReX_MLNodeLaplacian_sten.cpp
   MLMG/AMReX_MLNodeLaplacian_misc.cpp
   MLMG/AMReX_MLNodeLap_K.H
   MLMG/AMReX_MLNodeLap_${AMReX_SPACEDIM}D_K.H
   MLMG/AMReX_MLNodeTensorLaplacian.H
   MLMG/AMReX_MLNodeTensorLaplacian.cpp
   MLMG/AMReX_MLNodeTensorLap_K.H
   MLMG/AMReX_MLNodeTensorLap_${AMReX_SPACEDIM}D_K.H
   MLMG/AMReX_MLTensorOp.H
   MLMG/AMReX_MLTensorOp.cpp
   MLMG/AMReX_MLTensorOp_grad.cpp
   MLMG/AMReX_MLTensor_K.H
   MLMG/AMReX_MLTensor_${AMReX_SPACEDIM}D_K.H
   MLMG/AMReX_MLEBNodeFDLaplacian.H
   MLMG/AMReX_MLEBNodeFDLaplacian.cpp
   MLMG/AMReX_MLEBNodeFDLap_K.H
   MLMG/AMReX_MLEBNodeFDLap_${AMReX_SPACEDIM}D_K.H
   )

if (AMReX_SPACEDIM EQUAL 3)
   target_sources(amrex
      PRIVATE
      MLMG/AMReX_MLMG_2D_K.H
      MLMG/AMReX_MLPoisson_2D_K.H
      MLMG/AMReX_MLALap_2D_K.H
      )
endif ()

if (AMReX_EB)
   target_sources(amrex
      PRIVATE
      MLMG/AMReX_MLNodeLaplacian_eb.cpp
      MLMG/AMReX_MLEBABecLap.H
      MLMG/AMReX_MLEBABecLap.cpp
      MLMG/AMReX_MLEBABecLap_F.cpp
      MLMG/AMReX_MLEBABecLap_K.H
      MLMG/AMReX_MLEBABecLap_${AMReX_SPACEDIM}D_K.H
      MLMG/AMReX_MLEBTensorOp.H
      MLMG/AMReX_MLEBTensorOp.cpp
      MLMG/AMReX_MLEBTensorOp_bc.cpp
      MLMG/AMReX_MLEBTensor_K.H
      MLMG/AMReX_MLEBTensor_${AMReX_SPACEDIM}D_K.H
      )
endif ()

if (AMReX_FORTRAN)
   target_sources(amrex
      PRIVATE
      MLMG/AMReX_MLLinOp_F.H
      MLMG/AMReX_MLLinOp_nd.F90
      )
endif ()

if (AMReX_HYPRE)
   target_sources(amrex
      PRIVATE
      MLMG/AMReX_MLNodeLaplacian_hypre.cpp
      )
endif ()

if (AMReX_SPACEDIM EQUAL 3)

   target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/OpenBC>)

   target_sources(amrex
      PRIVATE
      OpenBC/AMReX_OpenBC.H
      OpenBC/AMReX_OpenBC_K.H
      OpenBC/AMReX_OpenBC.cpp
      )
endif ()
