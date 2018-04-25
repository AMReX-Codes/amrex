#
# Sources in folder MGML
# 
set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")


list ( APPEND ALLHEADERS AMReX_MLMG.H )
list ( APPEND CXXSRC     AMReX_MLMG.cpp )
list ( APPEND ALLHEADERS AMReX_MLMG_F.H )
list ( APPEND F90SRC     AMReX_MLMG_${DIM}d.F90 )

list ( APPEND ALLHEADERS AMReX_MLMGBndry.H )
list ( APPEND CXXSRC     AMReX_MLMGBndry.cpp )

list ( APPEND ALLHEADERS AMReX_MLLinOp.H )
list ( APPEND CXXSRC     AMReX_MLLinOp.cpp )
list ( APPEND ALLHEADERS AMReX_MLLinOp_F.H )
list ( APPEND F90SRC     AMReX_MLLinOp_nd.F90 )
list ( APPEND F90SRC     AMReX_MLLinOp_${DIM}d.F90 )

list ( APPEND ALLHEADERS AMReX_MLCellLinOp.H )
list ( APPEND CXXSRC     AMReX_MLCellLinOp.cpp )

list ( APPEND ALLHEADERS AMReX_MLNodeLinOp.H )
list ( APPEND CXXSRC     AMReX_MLNodeLinOp.cpp )

list ( APPEND ALLHEADERS AMReX_MLCGSolver.H )
list ( APPEND CXXSRC     AMReX_MLCGSolver.cpp )

list ( APPEND ALLHEADERS AMReX_MLABecLaplacian.H )
list ( APPEND CXXSRC     AMReX_MLABecLaplacian.cpp )
list ( APPEND ALLHEADERS AMReX_MLABecLap_F.H )
list ( APPEND F90SRC     AMReX_MLABecLap_${DIM}d.F90 )

list ( APPEND ALLHEADERS AMReX_MLALaplacian.H )
list ( APPEND CXXSRC     AMReX_MLALaplacian.cpp )
list ( APPEND ALLHEADERS AMReX_MLALap_F.H )
list ( APPEND F90SRC     AMReX_MLALap_${DIM}d.F90 )

list ( APPEND ALLHEADERS AMReX_MLPoisson.H )
list ( APPEND CXXSRC     AMReX_MLPoisson.cpp )
list ( APPEND ALLHEADERS AMReX_MLPoisson_F.H )
list ( APPEND F90SRC     AMReX_MLPoisson_${DIM}d.F90 )

list ( APPEND ALLHEADERS AMReX_MLNodeLaplacian.H )
list ( APPEND CXXSRC     AMReX_MLNodeLaplacian.cpp )
list ( APPEND ALLHEADERS AMReX_MLNodeLap_F.H )
list ( APPEND F90SRC     AMReX_MLNodeLap_${DIM}d.F90 )
list ( APPEND F90SRC     AMReX_MLNodeLap_nd.F90 )

# 
# Collect sources
# 
set ( ALLSRC   ${CXXSRC} ${F90SRC} )

#
# Add absolute path (needed by target_sources)
# 
prepend ( ALLSRC     ${CMAKE_CURRENT_LIST_DIR}/MLMG )
prepend ( ALLHEADERS ${CMAKE_CURRENT_LIST_DIR}/MLMG )

target_sources ( amrex PRIVATE ${ALLSRC} ${ALLHEADERS} )

# 
# Accumulate headers 
#
accumulate_headers (ALLHEADERS)   

#
# Add to include dirs
#
accumulate_include_dirs (${CMAKE_CURRENT_LIST_DIR}/MLMG)

#-------------------------------------------------------------------------------
#
# Sources in folder C_CellMG
# 
set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")


list ( APPEND ALLHEADERS AMReX_ABecLaplacian.H   AMReX_CGSolver.H   AMReX_LinOp.H   )
list ( APPEND CXXSRC     AMReX_ABecLaplacian.cpp AMReX_CGSolver.cpp AMReX_LinOp.cpp )

list ( APPEND ALLHEADERS AMReX_MultiGrid.H   AMReX_Laplacian.H )
list ( APPEND CXXSRC     AMReX_MultiGrid.cpp AMReX_Laplacian.cpp )

list ( APPEND ALLHEADERS AMReX_ABec_F.H AMReX_LO_F.H AMReX_LP_F.H AMReX_MG_F.H )
list ( APPEND F90SRC     AMReX_ABec_${DIM}D.F90 AMReX_ABec_UTIL.F90 AMReX_LO_${DIM}D.F90 )
list ( APPEND F90SRC     AMReX_LP_${DIM}D.F90 AMReX_MG_${DIM}D.F90 )

# 
# Collect sources
# 
set ( ALLSRC   ${CXXSRC} ${F90SRC} )

#
# Add absolute path (needed by target_sources)
# 
prepend ( ALLSRC     ${CMAKE_CURRENT_LIST_DIR}/C_CellMG )
prepend ( ALLHEADERS ${CMAKE_CURRENT_LIST_DIR}/C_CellMG )

target_sources ( amrex PRIVATE ${ALLSRC} ${ALLHEADERS} )

# 
# Accumulate headers 
#
accumulate_headers (ALLHEADERS)   

#
# Add to include dirs
#
accumulate_include_dirs (${CMAKE_CURRENT_LIST_DIR}/C_CellMG)


#-------------------------------------------------------------------------------
#
# Source in folder C_CellMG4
#
set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")

list ( APPEND ALLHEADERS AMReX_ABec2.H   AMReX_ABec4.H )
list ( APPEND CXXSRC     AMReX_ABec2.cpp AMReX_ABec4.cpp )

list ( APPEND ALLHEADERS AMReX_ABec2_F.H AMReX_ABec4_F.H )
list ( APPEND F90SRC     AMReX_ABec2_${DIM}D.F90 AMReX_ABec4_${DIM}D.F90 )


# 
# Collect sources
# 
set ( ALLSRC   ${CXXSRC} ${F90SRC} )

#
# Add absolute path (needed by target_sources)
# 
prepend ( ALLSRC     ${CMAKE_CURRENT_LIST_DIR}/C_CellMG4 )
prepend ( ALLHEADERS ${CMAKE_CURRENT_LIST_DIR}/C_CellMG4 )

target_sources ( amrex PRIVATE ${ALLSRC} ${ALLHEADERS} )

# 
# Accumulate headers 
#
accumulate_headers (ALLHEADERS)   

#
# Add to include dirs
#
accumulate_include_dirs ( ${CMAKE_CURRENT_LIST_DIR}/C_CellMG4 )

#
# Sources in folder C_TensorMG
# 
set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")


list ( APPEND ALLHEADERS AMReX_DivVis.H   AMReX_MCCGSolver.H   AMReX_MCInterpBndryData.H )
list ( APPEND CXXSRC     AMReX_DivVis.cpp AMReX_MCCGSolver.cpp AMReX_MCInterpBndryData.cpp )

list ( APPEND ALLHEADERS AMReX_MCLinOp.H   AMReX_MCMultiGrid.H )
list ( APPEND CXXSRC     AMReX_MCLinOp.cpp AMReX_MCMultiGrid.cpp )

list ( APPEND ALLHEADERS AMReX_DivVis_F.H AMReX_MCINTERPBNDRYDATA_F.H AMReX_MCLO_F.H )
list ( APPEND F77SRC     AMReX_MCLO_${DIM}D.F AMReX_MCINTERPBNDRYDATA_${DIM}D.F )

if (${DIM} EQUAL 2)
   list ( APPEND F77SRC AMReX_DV_2D.F )
elseif (${DIM} EQUAL 3)
   list ( APPEND F77SRC AMReX_DV_3D1.F AMReX_DV_3D2.F AMReX_DV_3D3.F )
else ()
   message ( AUTHOR_WARNING  "No specific source file for case DIM=1 in Src/LinearSolvers/C_TensorMG" )
endif ()

list ( APPEND F90SRC  amrex_tmg_util.F90 )

# 
# Collect sources
# 
set ( ALLSRC   ${CXXSRC} ${F77SRC} ${F90SRC} )

#
# Add absolute path (needed by target_sources)
# 
prepend ( ALLSRC     ${CMAKE_CURRENT_LIST_DIR}/C_TensorMG )
prepend ( ALLHEADERS ${CMAKE_CURRENT_LIST_DIR}/C_TensorMG )

target_sources ( amrex PRIVATE ${ALLSRC} ${ALLHEADERS} )

# 
# Accumulate headers 
#
accumulate_headers (ALLHEADERS)   

#
# Add to include dirs
#
accumulate_include_dirs ( ${CMAKE_CURRENT_LIST_DIR}/C_TensorMG)


#
# Add Fortran MG and C/Fortran interfaces ONLY IF FBaselib is enabled
# In the GNUMake system, there are two makefiles, F_MG/namely FParallelMG.mak
# and F_MG/GPackage.mak. Since I do not know exactly what they are doing, I am
# adding all the source files and skip the conditionals in the makefile
# ( and all the preprocessor directives enabled in those files ). 
# 
if (NOT ENABLE_FBASELIB)
   return ()
endif ()

#-------------------------------------------------------------------------------

#
# Source in folder F_MG
#
set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")



list ( APPEND CXXSRC  C_to_F_MG/AMReX_FMultiGrid.cpp  C_to_F_MG/AMReX_MGT_Solver.cpp )

list ( APPEND F90SRC
   F_MG/cc_applyop.f90  F_MG/cc_stencil.f90  F_MG/mg_tower.f90 F_MG/nodal_divu.f90
   F_MG/nodal_smoothers.f90  F_MG/stencil_types.f90  F_MG/cc_interface_stencil.f90
   F_MG/cc_stencil_fill.f90  F_MG/ml_cc.f90 F_MG/nodal_enforce_dirichlet_rhs.f90
   F_MG/nodal_stencil_apply.f90  F_MG/stencil_util.f90 F_MG/cc_mg_cpp.f90
   F_MG/coarsen_coeffs.f90  F_MG/ml_nd.f90 F_MG/nodal_interface_stencil.f90
   F_MG/nodal_stencil.f90  F_MG/tridiag.f90 F_MG/cc_mg_tower_smoother.f90
   F_MG/compute_defect.f90  F_MG/ml_norm.f90 F_MG/nodal_mask.f90
   F_MG/nodal_stencil_fill.f90 F_MG/cc_ml_resid.f90 F_MG/itsol.f90
   F_MG/ml_prolongation.f90  F_MG/nodal_mg_cpp.f90  F_MG/nodal_sum.f90 F_MG/cc_smoothers.f90
   F_MG/mg.f90 F_MG/ml_solve.f90 F_MG/nodal_mg_tower_smoother.f90 F_MG/nodal_sync_resid.f90
   F_MG/cc_stencil_apply.f90  F_MG/mg_prolongation.f90 F_MG/nodal_applyop.f90
   F_MG/nodal_newu.f90 )
# This one does not compile and it was no included in original
# CMakeLists.txt --> F_MG/sparse_solve.f90 )

list ( APPEND ALLHEADERS  C_to_F_MG/AMReX_FMultiGrid.H  C_to_F_MG/AMReX_MGT_Solver.H
   C_to_F_MG/AMReX_stencil_types.H F_MG/mg_cpp_f.h )


# Collect sources
set ( ALLSRC   ${CXXSRC} ${F90SRC}  )

#
# Add absolute path (needed by target_sources)
# 
prepend ( ALLSRC     ${CMAKE_CURRENT_LIST_DIR} )
prepend ( ALLHEADERS ${CMAKE_CURRENT_LIST_DIR} )

target_sources ( amrex PRIVATE ${ALLSRC} ${ALLHEADERS} )

# 
# Accumulate headers 
#
accumulate_headers (ALLHEADERS)   


#
# Add to include dirs
#
accumulate_include_dirs (${CMAKE_CURRENT_LIST_DIR}/C_to_F_MG)
accumulate_include_dirs (${CMAKE_CURRENT_LIST_DIR}/F_MG)
