#
# Sources in subdirectory MLMG
#
target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/MLMG>)

target_sources(amrex
   PRIVATE
   MLMG/AMReX_MLMG.H
   MLMG/AMReX_MLMG.cpp
   MLMG/AMReX_MLMG_F.H
   MLMG/AMReX_MLMG_${DIM}d.F90
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
   MLMG/AMReX_MLNodeLap_F.H
   MLMG/AMReX_MLNodeLap_${DIM}d.F90
   MLMG/AMReX_MLNodeLap_nd.F90
   MLMG/AMReX_MLTensorOp.H
   MLMG/AMReX_MLTensorOp.cpp
   MLMG/AMReX_MLTensor_K.H
   MLMG/AMReX_MLTensor_${DIM}D_K.H
   MLMG/AMReX_MacProjector.H
   MLMG/AMReX_MacProjector.cpp
   MLMG/AMReX_MG_K.H
   MLMG/AMReX_MG_${DIM}D_K.H
   )

if (ENABLE_EB)
   target_sources(amrex PRIVATE 
      MLMG/AMReX_MLEBABecLap.H
      MLMG/AMReX_MLEBABecLap.cpp
      MLMG/AMReX_MLEBABecLap_F.H
      MLMG/AMReX_MLEBABecLap_${DIM}d.F90
      MLMG/AMReX_MLEBABecLap_nd.F90
      )
endif ()


#
# Sources in folder C_CellMG
#
target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/C_CellMG>)

target_sources(amrex PRIVATE 
   C_CellMG/AMReX_ABecLaplacian.H
   C_CellMG/AMReX_ABecLaplacian.cpp
   C_CellMG/AMReX_CGSolver.H   
   C_CellMG/AMReX_CGSolver.cpp
   C_CellMG/AMReX_LinOp.H  
   C_CellMG/AMReX_LinOp.cpp
   C_CellMG/AMReX_MultiGrid.H   
   C_CellMG/AMReX_MultiGrid.cpp
   C_CellMG/AMReX_Laplacian.H
   C_CellMG/AMReX_Laplacian.cpp
   C_CellMG/AMReX_ABec_F.H
   C_CellMG/AMReX_LO_F.H
   C_CellMG/AMReX_LP_F.H
   C_CellMG/AMReX_MG_F.H
   C_CellMG/AMReX_ABec_${DIM}D.F90
   C_CellMG/AMReX_ABec_UTIL.F90
   C_CellMG/AMReX_LO_${DIM}D.F90
   C_CellMG/AMReX_LP_${DIM}D.F90
   C_CellMG/AMReX_MG_${DIM}D.F90
   )
