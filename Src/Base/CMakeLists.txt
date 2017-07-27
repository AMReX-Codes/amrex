set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")

set ( CXXSRC
   AMReX_Arena.cpp AMReX_BArena.cpp AMReX_BaseFab.cpp AMReX_BCRec.cpp
   AMReX_BLBackTrace.cpp AMReX_BoxArray.cpp AMReX_Box.cpp AMReX_BoxIterator.cpp
   AMReX_BoxDomain.cpp AMReX.cpp AMReX_BoxList.cpp AMReX_CArena.cpp
   AMReX_CoordSys.cpp AMReX_DistributionMapping.cpp AMReX_FabArrayBase.cpp
   AMReX_FabConv.cpp AMReX_FArrayBox.cpp AMReX_FPC.cpp AMReX_Geometry.cpp
   AMReX_MultiFabUtil.cpp AMReX_IArrayBox.cpp AMReX_IndexType.cpp
   AMReX_IntVect.cpp AMReX_iMultiFab.cpp AMReX_MemPool.cpp AMReX_MultiFab.cpp
   AMReX_MFCopyDescriptor.cpp AMReX_MFIter.cpp AMReX_NFiles.cpp
   AMReX_Orientation.cpp AMReX_ParallelDescriptor.cpp AMReX_ParmParse.cpp
   AMReX_Periodicity.cpp AMReX_PhysBCFunct.cpp AMReX_PlotFileUtil.cpp
   AMReX_RealBox.cpp AMReX_RealVect.cpp AMReX_Utility.cpp AMReX_VisMF.cpp
   AMReX_parmparse_fi.cpp AMReX_error_fi.cpp AMReX_BLProfiler.cpp
   AMReX_TinyProfiler.cpp )

set ( F77SRC
   AMReX_BLProfiler_F.f AMReX_BLBoxLib_F.f AMReX_bl_flush.f
   AMReX_BLParmParse_F.f AMReX_BLutil_F.f 
   AMReX_COORDSYS_${BL_SPACEDIM}D.F AMReX_FILCC_${BL_SPACEDIM}D.F )


set ( F90SRC
   AMReX_mempool_f.f90 AMReX_MultiFabUtil_${BL_SPACEDIM}d.f90
   AMReX_BaseFab_nd.f90 AMReX_filcc_mod.F90 AMReX_omp_mod.F90 AMReX_fort_mod.F90
   AMReX_parmparse_mod.F90 AMReX_string_mod.F90 AMReX_error_mod.F90 AMReX_bc_types_mod.F90 )

set ( ALLHEADERS
   AMReX_Arena.H AMReX_Array.H AMReX_ArrayLim.H AMReX_BArena.H AMReX_BaseFab.H
   AMReX_BCRec.H AMReX_BC_TYPES.H AMReX_BLassert.H AMReX_BLBackTrace.H AMReX_BLFort.H
   AMReX_BLProfiler.H AMReX_BoxArray.H AMReX_BoxDomain.H AMReX_Box.H AMReX_BoxIterator.H
   AMReX.H AMReX_BoxList.H AMReX_CArena.H AMReX_ccse-mpi.H AMReX_CONSTANTS.H
   AMReX_CoordSys.H AMReX_DistributionMapping.H AMReX_FabArray.H AMReX_FabArrayBase.H
   AMReX_FACopyDescriptor.H AMReX_FabConv.H AMReX_FabFactory.H AMReX_FArrayBox.H
   AMReX_FPC.H AMReX_Geometry.H AMReX_MultiFabUtil.H AMReX_IArrayBox.H
   AMReX_IndexType.H AMReX_IntVect.H AMReX_Looping.H AMReX_iMultiFab.H AMReX_MakeType.H
   AMReX_MemPool.H AMReX_MultiFab.H AMReX_MFCopyDescriptor.H AMReX_MFIter.H AMReX_NFiles.H
   AMReX_Orientation.H AMReX_ParallelDescriptor.H AMReX_ParmParse.H AMReX_Array.H
   AMReX_Periodicity.H AMReX_PlotFileUtil.H AMReX_Print.H AMReX_RealBox.H AMReX_RealVect.H
   AMReX_REAL.H AMReX_SPACE.H AMReX_Tuple.H AMReX_TypeTraits.H AMReX_Utility.H
   AMReX_parstream.H AMReX_VisMF.H AMReX_PhysBCFunct.H AMReX_SPMD.H AMReX_SPMDI.H
   AMReX_TinyProfiler.H AMReX_BCUtil.H
   # Fortran Headers
   AMReX_bc_types.fi AMReX_COORDSYS_F.H AMReX_SPACE_F.H AMReX_BaseFab_f.H AMReX_filcc_f.H )



# Accumulate sources
set ( ALLSRC ${CXXSRC} ${F90SRC} ${F77SRC} )

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
