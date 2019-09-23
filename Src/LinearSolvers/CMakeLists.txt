#
# Sources in subdirectory MLMG
#
target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/MLMG>)

target_sources(amrex
   PRIVATE
   MLMG/AMReX_MLMG.H
   MLMG/AMReX_MLMG.cpp
   MLMG/AMReX_MLMG_K.H
   MLMG/AMReX_MLMG_${DIM}D_K.H
   MLMG/AMReX_MLMGBndry.H
   MLMG/AMReX_MLMGBndry.cpp
   MLMG/AMReX_MLLinOp.H
   MLMG/AMReX_MLLinOp.cpp
   MLMG/AMReX_MLLinOp_K.H
   MLMG/AMReX_MLLinOp_F.H
   MLMG/AMReX_MLLinOp_nd.F90
   MLMG/AMReX_MLCellLinOp.H
   MLMG/AMReX_MLCellLinOp.cpp
   MLMG/AMReX_MLNodeLinOp.H
   MLMG/AMReX_MLNodeLinOp.cpp
   MLMG/AMReX_MLCellABecLap.H
   MLMG/AMReX_MLCellABecLap.cpp
   MLMG/AMReX_MLCGSolver.H
   MLMG/AMReX_MLCGSolver.cpp
   MLMG/AMReX_MLABecLaplacian.H
   MLMG/AMReX_MLABecLaplacian.cpp
   MLMG/AMReX_MLABecLap_K.H
   MLMG/AMReX_MLABecLap_${DIM}D_K.H
   MLMG/AMReX_MLALaplacian.H
   MLMG/AMReX_MLALaplacian.cpp
   MLMG/AMReX_MLALap_K.H
   MLMG/AMReX_MLALap_${DIM}D_K.H
   MLMG/AMReX_MLPoisson.H
   MLMG/AMReX_MLPoisson.cpp
   MLMG/AMReX_MLPoisson_K.H
   MLMG/AMReX_MLPoisson_${DIM}D_K.H
   MLMG/AMReX_MLNodeLaplacian.H
   MLMG/AMReX_MLNodeLaplacian.cpp
   MLMG/AMReX_MLNodeLap_K.H
   MLMG/AMReX_MLNodeLap_${DIM}D_K.H
   MLMG/AMReX_MLNodeLap_F.H
   MLMG/AMReX_MLNodeLap_${DIM}d.F90
   MLMG/AMReX_MLTensorOp.H
   MLMG/AMReX_MLTensorOp.cpp
   MLMG/AMReX_MLTensor_K.H
   MLMG/AMReX_MLTensor_${DIM}D_K.H
   MLMG/AMReX_MacProjector.H
   MLMG/AMReX_MacProjector.cpp
   )

if (ENABLE_EB)
   target_sources(amrex PRIVATE 
      MLMG/AMReX_MLEBABecLap.H
      MLMG/AMReX_MLEBABecLap.cpp
      MLMG/AMReX_MLEBABecLap_K.H
      MLMG/AMReX_MLEBABecLap_${DIM}D_K.H
      MLMG/AMReX_MLEBTensorOp.H
      MLMG/AMReX_MLEBTensorOp.cpp
      MLMG/AMReX_MLEBTensor_K.H
      MLMG/AMReX_MLEBTensor_${DIM}D_K.H
      )
endif ()
