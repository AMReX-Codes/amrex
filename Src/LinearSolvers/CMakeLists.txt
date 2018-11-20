add_sources ( MLMG/AMReX_MLMG.H )
add_sources ( MLMG/AMReX_MLMG.cpp )
add_sources ( MLMG/AMReX_MLMG_F.H )
add_sources ( MLMG/AMReX_MLMG_${DIM}d.F90 )

add_sources ( MLMG/AMReX_MLMGBndry.H )
add_sources ( MLMG/AMReX_MLMGBndry.cpp )

add_sources ( MLMG/AMReX_MLLinOp.H )
add_sources ( MLMG/AMReX_MLLinOp.cpp )
add_sources ( MLMG/AMReX_MLLinOp_F.H )
add_sources ( MLMG/AMReX_MLLinOp_nd.F90 )
add_sources ( MLMG/AMReX_MLLinOp_${DIM}d.F90 )

add_sources ( MLMG/AMReX_MLCellLinOp.H )
add_sources ( MLMG/AMReX_MLCellLinOp.cpp )

add_sources ( MLMG/AMReX_MLNodeLinOp.H )
add_sources ( MLMG/AMReX_MLNodeLinOp.cpp )

add_sources ( MLMG/AMReX_MLCellABecLap.H )
add_sources ( MLMG/AMReX_MLCellABecLap.cpp )

add_sources ( MLMG/AMReX_MLCGSolver.H )
add_sources ( MLMG/AMReX_MLCGSolver.cpp )

add_sources ( MLMG/AMReX_MLABecLaplacian.H )
add_sources ( MLMG/AMReX_MLABecLaplacian.cpp )
add_sources ( MLMG/AMReX_MLABecLap_F.H )
add_sources ( MLMG/AMReX_MLABecLap_${DIM}d.F90 )

if (ENABLE_EB)
    add_sources ( MLMG/AMReX_MLEBABecLap.H )
    add_sources ( MLMG/AMReX_MLEBABecLap.cpp )
    add_sources ( MLMG/AMReX_MLEBABecLap_F.H )
    add_sources ( MLMG/AMReX_MLEBABecLap_${DIM}d.F90 )
    add_sources ( MLMG/AMReX_MLEBABecLap_nd.F90 )
endif ()

add_sources ( MLMG/AMReX_MLALaplacian.H )
add_sources ( MLMG/AMReX_MLALaplacian.cpp )
add_sources ( MLMG/AMReX_MLALap_F.H )
add_sources ( MLMG/AMReX_MLALap_${DIM}d.F90 )

add_sources ( MLMG/AMReX_MLPoisson.H )
add_sources ( MLMG/AMReX_MLPoisson.cpp )
add_sources ( MLMG/AMReX_MLPoisson_F.H )
add_sources ( MLMG/AMReX_MLPoisson_${DIM}d.F90 )

add_sources ( MLMG/AMReX_MLNodeLaplacian.H )
add_sources ( MLMG/AMReX_MLNodeLaplacian.cpp )
add_sources ( MLMG/AMReX_MLNodeLap_F.H )
add_sources ( MLMG/AMReX_MLNodeLap_${DIM}d.F90 )
add_sources ( MLMG/AMReX_MLNodeLap_nd.F90 )

add_sources ( MLMG/AMReX_MacProjector.H )
add_sources ( MLMG/AMReX_MacProjector.cpp )


#-------------------------------------------------------------------------------
#
# Sources in folder C_CellMG
# 
add_sources ( C_CellMG/AMReX_ABecLaplacian.H   C_CellMG/AMReX_CGSolver.H   C_CellMG/AMReX_LinOp.H   )
add_sources ( C_CellMG/AMReX_ABecLaplacian.cpp C_CellMG/AMReX_CGSolver.cpp C_CellMG/AMReX_LinOp.cpp )

add_sources ( C_CellMG/AMReX_MultiGrid.H   C_CellMG/AMReX_Laplacian.H )
add_Sources ( C_CellMG/AMReX_MultiGrid.cpp C_CellMG/AMReX_Laplacian.cpp )

add_sources ( C_CellMG/AMReX_ABec_F.H C_CellMG/AMReX_LO_F.H C_CellMG/AMReX_LP_F.H C_CellMG/AMReX_MG_F.H )
add_sources ( C_CellMG/AMReX_ABec_${DIM}D.F90 C_CellMG/AMReX_ABec_UTIL.F90 C_CellMG/AMReX_LO_${DIM}D.F90 )
add_sources ( C_CellMG/AMReX_LP_${DIM}D.F90 C_CellMG/AMReX_MG_${DIM}D.F90 )

