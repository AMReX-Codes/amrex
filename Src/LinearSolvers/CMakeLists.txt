set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")


set (CXXSRC
   C_CellMG/AMReX_ABecLaplacian.cpp  C_CellMG/AMReX_CGSolver.cpp  C_CellMG/AMReX_Laplacian.cpp
   C_CellMG/AMReX_LinOp.cpp  C_CellMG/AMReX_MultiGrid.cpp
   C_CellMG4/AMReX_ABec2.cpp  C_CellMG4/AMReX_ABec4.cpp
   C_TensorMG/AMReX_DivVis.cpp  C_TensorMG/AMReX_MCCGSolver.cpp
   C_TensorMG/AMReX_MCInterpBndryData.cpp  C_TensorMG/AMReX_MCLinOp.cpp
   C_TensorMG/AMReX_MCMultiGrid.cpp )

set (F77SRC
   C_CellMG/AMReX_ABec_${BL_SPACEDIM}D.F  C_CellMG/AMReX_LO_${BL_SPACEDIM}D.F
   C_CellMG/AMReX_LP_${BL_SPACEDIM}D.F  C_CellMG/AMReX_MG_${BL_SPACEDIM}D.F 
   C_CellMG/AMReX_ABec_UTIL.F
   C_CellMG4/AMReX_ABec2_${BL_SPACEDIM}D.F  C_CellMG4/AMReX_ABec4_${BL_SPACEDIM}D.F
   C_TensorMG/AMReX_MCINTERPBNDRYDATA_${BL_SPACEDIM}D.F
   )

if ( BL_SPACEDIM EQUAL 2)
   list (APPEND F77SRC C_TensorMG/AMReX_DV_2D.F)
elseif (BL_SPACEDIM EQUAL 3)
   list (APPEND F77SRC C_TensorMG/AMReX_DV_3D1.F C_TensorMG/AMReX_DV_3D2.F
      C_TensorMG/AMReX_DV_3D3.F)
endif ()
   
   
set (ALLHEADERS
   C_CellMG/AMReX_ABec_F.H  C_CellMG/AMReX_CGSolver.H   C_CellMG/AMReX_LinOp.H
   C_CellMG/AMReX_LP_F.H    C_CellMG/AMReX_MultiGrid.H  C_CellMG/AMReX_ABecLaplacian.H
   C_CellMG/AMReX_Laplacian.H  C_CellMG/AMReX_LO_F.H    C_CellMG/AMReX_MG_F.H
   C_CellMG4/AMReX_ABec2_F.H  C_CellMG4/AMReX_ABec2.H  C_CellMG4/AMReX_ABec4_F.H
   C_CellMG4/AMReX_ABec4.H
   C_TensorMG/AMReX_DivVis_F.H  C_TensorMG/AMReX_MCCGSolver.H
   C_TensorMG/AMReX_MCInterpBndryData.H  C_TensorMG/AMReX_MCLO_F.H
   C_TensorMG/AMReX_DivVis.H    C_TensorMG/AMReX_MCINTERPBNDRYDATA_F.H
   C_TensorMG/AMReX_MCLinOp.H   C_TensorMG/AMReX_MCMultiGrid.H
   C_to_F_MG/AMReX_FMultiGrid.H  C_to_F_MG/AMReX_MGT_Solver.H  C_to_F_MG/AMReX_stencil_types.H
   F_MG/mg_cpp_f.h )


#
# Add Fortran MG and C/Fortran interfaces is FBaselib is enabled
# 
if (ENABLE_FBASELIB)

   list ( APPEND CXXSRC  C_to_F_MG/AMReX_FMultiGrid.cpp  C_to_F_MG/AMReX_MGT_Solver.cpp )

   set (F90SRC
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
   
endif ()

# Collect sources
set ( ALLSRC   ${CXXSRC} ${F90SRC} ${F77SRC} )

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

