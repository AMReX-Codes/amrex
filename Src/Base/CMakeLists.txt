set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")


#
# Utility classes
#
list ( APPEND ALLHEADERS AMReX_ccse-mpi.H AMReX_Array.H AMReX_Vector.H AMReX_BaseIndex.H AMReX_Tuple.H )

list ( APPEND ALLHEADERS AMReX.H)
list ( APPEND CXXSRC     AMReX.cpp AMReX_error_fi.cpp )

list ( APPEND CXXSRC     AMReX_ParmParse.cpp AMReX_parmparse_fi.cpp AMReX_Utility.cpp )
list ( APPEND ALLHEADERS AMReX_ParmParse.H AMReX_Utility.H AMReX_BLassert.H AMReX_ArrayLim.H )

list ( APPEND ALLHEADERS AMReX_REAL.H AMReX_CONSTANTS.H AMReX_SPACE.H AMReX_SPACE_F.H )

list ( APPEND CXXSRC     AMReX_DistributionMapping.cpp AMReX_ParallelDescriptor.cpp )
list ( APPEND ALLHEADERS AMReX_DistributionMapping.H AMReX_ParallelDescriptor.H )

list ( APPEND ALLHEADERS AMReX_ParallelReduce.H )

list ( APPEND CXXSRC     AMReX_VisMF.cpp AMReX_Arena.cpp AMReX_BArena.cpp AMReX_CArena.cpp )
list ( APPEND ALLHEADERS AMReX_VisMF.H AMReX_Arena.H AMReX_BArena.H AMReX_CArena.H )

list ( APPEND ALLHEADERS AMReX_BLProfiler.H AMReX_BLBackTrace.H AMReX_BLFort.H )

list ( APPEND CXXSRC     AMReX_NFiles.cpp )
list ( APPEND ALLHEADERS AMReX_NFiles.H )
   
list ( APPEND CXXSRC     AMReX_SPMD.cpp )
list ( APPEND ALLHEADERS AMReX_SPMD.H AMReX_SPMDI.H )

list ( APPEND CXXSRC     AMReX_parstream.cpp )
list ( APPEND ALLHEADERS AMReX_parstream.H )

#
# FAB I/O stuff
# 
list ( APPEND CXXSRC     AMReX_FabConv.cpp AMReX_FPC.cpp )
list ( APPEND ALLHEADERS AMReX_FabConv.H AMReX_FPC.H AMReX_Print.H )

#
# Index space
#
list ( APPEND CXXSRC     AMReX_Box.cpp AMReX_BoxIterator.cpp AMReX_IntVect.cpp
   AMReX_IndexType.cpp AMReX_Orientation.cpp AMReX_Periodicity.cpp )
list ( APPEND ALLHEADERS AMReX_Box.H AMReX_BoxIterator.H AMReX_IntVect.H
   AMReX_IndexType.H AMReX_Orientation.H AMReX_Periodicity.H  )

#
# Real space
#
list ( APPEND CXXSRC     AMReX_RealBox.cpp AMReX_RealVect.cpp )
list ( APPEND ALLHEADERS AMReX_RealBox.H   AMReX_RealVect.H )

#
# Unions of rectangle
#
list ( APPEND CXXSRC     AMReX_BoxList.cpp AMReX_BoxArray.cpp AMReX_BoxDomain.cpp )
list ( APPEND ALLHEADERS AMReX_BoxList.H AMReX_BoxArray.H AMReX_BoxDomain.H )

#
# Fortran array data
# 
list ( APPEND CXXSRC     AMReX_FArrayBox.cpp AMReX_IArrayBox.cpp AMReX_BaseFab.cpp )
list ( APPEND ALLHEADERS AMReX_FArrayBox.H AMReX_IArrayBox.H AMReX_Looping.H AMReX_MakeType.H
   AMReX_TypeTraits.H AMReX_BaseFab.H AMReX_FabFactory.H )

#
# Fortran data defined on unions of rectangles.
#
list ( APPEND CXXSRC     AMReX_MultiFab.cpp AMReX_MFCopyDescriptor.cpp )
list ( APPEND ALLHEADERS AMReX_MultiFab.H AMReX_MFCopyDescriptor.H )

list ( APPEND CXXSRC     AMReX_iMultiFab.cpp )
list ( APPEND ALLHEADERS AMReX_iMultiFab.H )

list ( APPEND CXXSRC     AMReX_FabArrayBase.cpp AMReX_MFIter.cpp )
list ( APPEND ALLHEADERS AMReX_FabArray.H AMReX_FACopyDescriptor.H )
list ( APPEND ALLHEADERS AMReX_FabArrayBase.H AMReX_MFIter.H AMReX_LayoutData.H)

#
# Geometry / Coordinate system routines.
# In GNUMake system, this is included only if BL_NO_FORT=FALSE 
# 
list ( APPEND CXXSRC     AMReX_CoordSys.cpp AMReX_Geometry.cpp AMReX_MultiFabUtil.cpp )
list ( APPEND ALLHEADERS AMReX_CoordSys.H AMReX_Geometry.H AMReX_MultiFabUtil.H )

#
# Boundary-related
#
list ( APPEND CXXSRC     AMReX_BCRec.cpp AMReX_PhysBCFunct.cpp AMReX_BCUtil.cpp )
list ( APPEND ALLHEADERS AMReX_BCRec.H AMReX_PhysBCFunct.H AMReX_BC_TYPES.H AMReX_BCUtil.H )

#
# Plotfile
# 
list ( APPEND CXXSRC     AMReX_PlotFileUtil.cpp )
list ( APPEND ALLHEADERS AMReX_PlotFileUtil.H )

#
# Fortran interface routines.
# In GNUMake system, this is included only if BL_NO_FORT=FALSE
#
list ( APPEND F77SRC     AMReX_COORDSYS_${DIM}D.F )
list ( APPEND ALLHEADERS AMReX_COORDSYS_F.H AMReX_filcc_f.H )

list ( APPEND F77SRC     AMReX_BLutil_F.f AMReX_BLParmParse_F.f AMReX_BLBoxLib_F.f
   AMReX_BLProfiler_F.f)

list ( APPEND F90SRC     AMReX_BaseFab_nd.f90 AMReX_MultiFabUtil_${DIM}d.f90 AMReX_MultiFabUtil_nd.F90 AMReX_FILCC_${DIM}D.F90 )
list ( APPEND ALLHEADERS AMReX_BaseFab_f.H AMReX_MultiFabUtil_F.H )

list ( APPEND F90SRC     AMReX_filcc_mod.F90 AMReX_omp_mod.F90 AMReX_fort_mod.F90 AMReX_constants_mod.f90 )
list ( APPEND F90SRC     AMReX_error_mod.F90 AMReX_parmparse_mod.F90 AMReX_string_mod.F90 )
list ( APPEND F90SRC     AMReX_bc_types_mod.F90 )

list ( APPEND CXXSRC     AMReX_BLProfiler.cpp AMReX_BLBackTrace.cpp )


# LAZY mode not supported in CMake yet
# ifeq ($(LAZY),TRUE)
# C$(AMREX_BASE)_sources += AMReX_Lazy.cpp
# C$(AMREX_BASE)_headers += AMReX_Lazy.H
# endif

#
# Memory pool
#
list ( APPEND CXXSRC     AMReX_MemPool.cpp )
list ( APPEND ALLHEADERS AMReX_MemPool.H )
list ( APPEND F90SRC     AMReX_mempool_f.f90 ) # if BL_NO_FORT = FALSE

# UPCXX mode not supperted in CMake yet
# ifeq ($(USE_UPCXX),TRUE)
#   C$(AMREX_BASE)_sources += AMReX_BLPgas.cpp
#  C$(AMREX_BASE)_headers += AMReX_BLPgas.H
# endif

#
# Memory Profiler
#
if (ENABLE_MEM_PROFILE)
   list ( APPEND CXXSRC     AMReX_MemProfiler.cpp )
   list ( APPEND ALLHEADERS AMReX_MemProfiler.H )
endif ()

#
# Tiny Profiler
#
if (ENABLE_TINY_PROFILE)
   list ( APPEND CXXSRC     AMReX_TinyProfiler.cpp )
   list ( APPEND ALLHEADERS AMReX_TinyProfiler.H )
endif ()

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

#
# Add to include dirs
#
accumulate_include_dirs (${CMAKE_CURRENT_LIST_DIR})
