#
# Utility classes
#
add_sources( AMReX_ccse-mpi.H AMReX_Array.H AMReX_Vector.H AMReX_BaseIndex.H AMReX_Tuple.H )
add_sources( AMReX.cpp AMReX_error_fi.cpp )

add_sources( AMReX.H AMReX_Exception.H )
add_sources( AMReX.cpp AMReX_error_fi.cpp )

add_sources( AMReX_Extension.H )

add_sources( AMReX_IndexSequence.H )

add_sources( AMReX_ParmParse.cpp AMReX_parmparse_fi.cpp AMReX_Utility.cpp )
add_sources( AMReX_ParmParse.H AMReX_Utility.H AMReX_BLassert.H AMReX_ArrayLim.H )

add_sources( AMReX_REAL.H AMReX_CONSTANTS.H AMReX_SPACE.H )

add_sources( AMReX_DistributionMapping.cpp AMReX_ParallelDescriptor.cpp )
add_sources( AMReX_DistributionMapping.H AMReX_ParallelDescriptor.H )

add_sources( AMReX_ParallelReduce.H )

add_sources( AMReX_ForkJoin.H AMReX_ParallelContext.H )
add_sources( AMReX_ForkJoin.cpp AMReX_ParallelContext.cpp )

add_sources( AMReX_VisMF.cpp AMReX_Arena.cpp AMReX_BArena.cpp AMReX_CArena.cpp )
add_sources( AMReX_VisMF.H AMReX_Arena.H AMReX_BArena.H AMReX_CArena.H )

add_sources( AMReX_BLProfiler.H AMReX_BLBackTrace.H AMReX_BLFort.H )

add_sources( AMReX_NFiles.cpp )
add_sources( AMReX_NFiles.H )
   
add_sources( AMReX_parstream.cpp )
add_sources( AMReX_parstream.H )

#
# I/O stuff
# 
add_sources( AMReX_FabConv.cpp AMReX_FPC.cpp AMReX_IntConv.cpp AMReX_VectorIO.cpp)
add_sources( AMReX_FabConv.H AMReX_FPC.H AMReX_Print.H AMReX_IntConv.H AMReX_VectorIO.H)

#
# Index space
#
add_sources( AMReX_Box.cpp AMReX_BoxIterator.cpp AMReX_IntVect.cpp
   AMReX_IndexType.cpp AMReX_Orientation.cpp AMReX_Periodicity.cpp )
add_sources( AMReX_Box.H AMReX_BoxIterator.H AMReX_IntVect.H
   AMReX_IndexType.H AMReX_Orientation.H AMReX_Periodicity.H  )

#
# Real space
#
add_sources( AMReX_RealBox.cpp AMReX_RealVect.cpp )
add_sources( AMReX_RealBox.H   AMReX_RealVect.H )

#
# Unions of rectangle
#
add_sources( AMReX_BoxList.cpp AMReX_BoxArray.cpp AMReX_BoxDomain.cpp )
add_sources( AMReX_BoxList.H AMReX_BoxArray.H AMReX_BoxDomain.H )

#
# Fortran array data
# 
add_sources( AMReX_FArrayBox.cpp AMReX_IArrayBox.cpp AMReX_BaseFab.cpp )
add_sources( AMReX_FArrayBox.H AMReX_IArrayBox.H AMReX_MakeType.H
   AMReX_TypeTraits.H AMReX_BaseFab.H AMReX_FabFactory.H )
add_sources( AMReX_BaseFabUtility.H )

#
# Fortran data defined on unions of rectangles.
#
add_sources( AMReX_MultiFab.cpp AMReX_MFCopyDescriptor.cpp )
add_sources( AMReX_MultiFab.H AMReX_MFCopyDescriptor.H )

add_sources( AMReX_iMultiFab.cpp )
add_sources( AMReX_iMultiFab.H )

add_sources( AMReX_FabArrayBase.cpp AMReX_MFIter.cpp )
add_sources( AMReX_FabArray.H AMReX_FACopyDescriptor.H AMReX_FabArrayCommI.H )
add_sources( AMReX_FabArrayUtility.H )
add_sources( AMReX_FabArrayBase.H AMReX_MFIter.H AMReX_LayoutData.H)

#
# Geometry / Coordinate system routines.
# In GNUMake system, this is included only if BL_NO_FORT=FALSE 
# 
add_sources( AMReX_CoordSys.cpp AMReX_Geometry.cpp AMReX_MultiFabUtil.cpp )
add_sources( AMReX_CoordSys.H AMReX_Geometry.H AMReX_MultiFabUtil.H )

#
# Boundary-related
#
add_sources( AMReX_BCRec.cpp AMReX_PhysBCFunct.cpp AMReX_BCUtil.cpp )
add_sources( AMReX_BCRec.H AMReX_PhysBCFunct.H AMReX_BC_TYPES.H AMReX_BCUtil.H )

#
# Plotfile
# 
add_sources( AMReX_PlotFileUtil.cpp )
add_sources( AMReX_PlotFileUtil.H )

#
# Fortran interface routines.
# In GNUMake system, this is included only if BL_NO_FORT=FALSE
#
add_sources( AMReX_COORDSYS_${DIM}D_C.H )
add_sources( AMReX_COORDSYS_C.H AMReX_filcc_f.H )

add_sources( AMReX_BLutil_F.F90 AMReX_BLBoxLib_F.F90
   AMReX_BLProfiler_F.F90)

add_sources( AMReX_FILCC_${DIM}D.F90 )
add_sources( AMReX_MultiFabUtil_${DIM}D_C.H AMReX_MultiFabUtil_nd_C.H )
add_sources( AMReX_MultiFabUtil_C.H )

add_sources( AMReX_FilCC_${DIM}D_C.H )
add_sources( AMReX_FilCC_C.H )

add_sources( AMReX_filcc_mod.F90 AMReX_omp_mod.F90 AMReX_acc_mod.F90 AMReX_fort_mod.F90 AMReX_constants_mod.f90 )
add_sources( AMReX_error_mod.F90 AMReX_parmparse_mod.F90 AMReX_string_mod.F90 )
add_sources( AMReX_bc_types_mod.F90 )
add_sources( AMReX_ParallelDescriptor_F.F90 )

add_sources( AMReX_io_mod.F90 )

add_sources( AMReX_BLProfiler.cpp AMReX_BLBackTrace.cpp )


# LAZY mode not supported in CMake yet
# ifeq ($(LAZY),TRUE)
# C$(AMREX_BASE)_sources += AMReX_Lazy.cpp
# C$(AMREX_BASE)_headers += AMReX_Lazy.H
# endif

#
# Memory pool
#
add_sources( AMReX_MemPool.cpp )
add_sources( AMReX_MemPool.H )
add_sources( AMReX_mempool_mod.F90 ) # if BL_NO_FORT = FALSE

#
# Memory Profiler
#
if (ENABLE_MEM_PROFILE)
   add_sources( AMReX_MemProfiler.cpp )
   add_sources( AMReX_MemProfiler.H )
endif ()

#
# Tiny Profiler
#
if (ENABLE_TINY_PROFILE)
   add_sources( AMReX_TinyProfiler.cpp )
   add_sources( AMReX_TinyProfiler.H )
endif ()

add_sources( AMReX_Gpu.H AMReX_GpuQualifiers.H AMReX_GpuControl.H AMReX_GpuLaunch.H AMReX_GpuError.H )
add_sources( AMReX_GpuControl.cpp )

#
# This is temporaty CUDA stuff
# For now we include the following files, but they should really just
# be included when ENABLE_CUDA is set to yes 
#
#if (ENABLE_CUDA)
add_sources( AMReX_CudaAllocators.H AMReX_CudaRange.H )
add_sources( AMReX_CudaMemory.H AMReX_CudaAllocators.H AMReX_CudaLaunch.H AMReX_CudaContainers.H)
add_sources( AMReX_CudaReduce.H )
add_sources( AMReX_CudaLaunch.cpp AMReX_CudaAllocators.cpp)
add_sources( AMReX_CudaDevice.H AMReX_CudaDevice.cpp )
add_sources( AMReX_CudaUtility.H AMReX_CudaUtility.cpp )
add_sources( AMReX_CudaAsyncFab.H AMReX_CudaAsyncFab.cpp AMReX_CudaAsyncFabImpl.H AMReX_CudaAsyncFabImpl.cpp )
add_sources( AMReX_CudaAsyncArray.H AMReX_CudaAsyncArray.cpp )
#endif()

