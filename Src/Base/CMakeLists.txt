target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( amrex
   PRIVATE
   # Utility classes ---------------------------------------------------------
   AMReX_ccse-mpi.H
   AMReX_Math.H
   AMReX_Algorithm.H
   AMReX_Array.H
   AMReX_BlockMutex.H
   AMReX_BlockMutex.cpp
   AMReX_GpuComplex.H
   AMReX_Vector.H
   AMReX_Tuple.H
   AMReX.cpp
   AMReX_error_fi.cpp
   AMReX.H
   AMReX_Exception.H
   AMReX.cpp
   AMReX_error_fi.cpp
   AMReX_Extension.H
   AMReX_IndexSequence.H
   AMReX_PODVector.H
   AMReX_ParmParse.cpp
   AMReX_parmparse_fi.cpp
   AMReX_ParmParse.H
   AMReX_Functional.H
   AMReX_Utility.H
   AMReX_Utility.cpp
   AMReX_FileSystem.H
   AMReX_FileSystem.cpp
   AMReX_Reduce.H
   AMReX_Scan.H
   AMReX_Partition.H
   AMReX_Random.H
   AMReX_Random.cpp
   AMReX_BLassert.H
   AMReX_ArrayLim.H
   AMReX_REAL.H
   AMReX_INT.H
   AMReX_CONSTANTS.H
   AMReX_SPACE.H
   AMReX_DistributionMapping.H
   AMReX_DistributionMapping.cpp
   AMReX_ParallelDescriptor.H
   AMReX_ParallelDescriptor.cpp
   AMReX_OpenMP.H
   AMReX_ParallelReduce.H
   AMReX_ForkJoin.H
   AMReX_ForkJoin.cpp
   AMReX_ParallelContext.H
   AMReX_ParallelContext.cpp
   AMReX_VisMF.H
   AMReX_VisMF.cpp
   AMReX_AsyncOut.H
   AMReX_AsyncOut.cpp
   AMReX_BackgroundThread.H
   AMReX_BackgroundThread.cpp
   AMReX_Arena.H
   AMReX_Arena.cpp
   AMReX_BArena.H
   AMReX_BArena.cpp
   AMReX_CArena.H
   AMReX_CArena.cpp
   AMReX_DArena.H
   AMReX_DArena.cpp
   AMReX_EArena.H
   AMReX_EArena.cpp
   AMReX_BLProfiler.H
   AMReX_BLBackTrace.H
   AMReX_BLFort.H
   AMReX_NFiles.H
   AMReX_NFiles.cpp
   AMReX_parstream.H
   AMReX_parstream.cpp
   # I/O stuff  --------------------------------------------------------------
   AMReX_FabConv.H
   AMReX_FabConv.cpp
   AMReX_FPC.H
   AMReX_FPC.cpp
   AMReX_VectorIO.H
   AMReX_VectorIO.cpp
   AMReX_Print.H
   AMReX_IntConv.H
   AMReX_IntConv.cpp
   # Index space -------------------------------------------------------------
   AMReX_Box.H
   AMReX_Box.cpp
   AMReX_BoxIterator.H
   AMReX_BoxIterator.cpp
   AMReX_Dim3.H
   AMReX_IntVect.H
   AMReX_IntVect.cpp
   AMReX_IndexType.H
   AMReX_IndexType.cpp
   AMReX_Loop.H
   AMReX_Orientation.H
   AMReX_Orientation.cpp
   AMReX_Periodicity.H
   AMReX_Periodicity.cpp
   # Real space --------------------------------------------------------------
   AMReX_RealBox.H
   AMReX_RealBox.cpp
   AMReX_RealVect.H
   AMReX_RealVect.cpp
   # Unions of rectangle -----------------------------------------------------
   AMReX_BoxList.H
   AMReX_BoxList.cpp
   AMReX_BoxArray.H
   AMReX_BoxArray.cpp
   AMReX_BoxDomain.H
   AMReX_BoxDomain.cpp
   # Fortran array data ------------------------------------------------------
   AMReX_FArrayBox.H
   AMReX_FArrayBox.cpp
   AMReX_IArrayBox.H
   AMReX_IArrayBox.cpp
   AMReX_BaseFab.H
   AMReX_BaseFab.cpp
   AMReX_Array4.H
   AMReX_MakeType.H
   AMReX_TypeTraits.H
   AMReX_FabFactory.H
   AMReX_BaseFabUtility.H
   # Fortran data defined on unions of rectangles ----------------------------
   AMReX_MultiFab.cpp
   AMReX_MultiFab.H
   AMReX_MFCopyDescriptor.cpp
   AMReX_MFCopyDescriptor.H
   AMReX_iMultiFab.cpp
   AMReX_iMultiFab.H
   AMReX_FabArrayBase.cpp
   AMReX_FabArrayBase.H
   AMReX_MFIter.cpp
   AMReX_MFIter.H
   AMReX_FabArray.H
   AMReX_FACopyDescriptor.H
   AMReX_FabArrayCommI.H
   AMReX_FBI.H
   AMReX_PCI.H
   AMReX_FabArrayUtility.H
   AMReX_LayoutData.H
   # Geometry / Coordinate system routines -----------------------------------
   AMReX_CoordSys.cpp
   AMReX_CoordSys.H
   AMReX_COORDSYS_${DIM}D_C.H
   AMReX_COORDSYS_C.H
   AMReX_Geometry.cpp
   AMReX_Geometry.H
   AMReX_MultiFabUtil.cpp
   AMReX_MultiFabUtil.H
   AMReX_MultiFabUtilI.H
   AMReX_MultiFabUtil_${DIM}D_C.H
   AMReX_MultiFabUtil_nd_C.H
   AMReX_MultiFabUtil_C.H
   # Boundary-related --------------------------------------------------------
   AMReX_BCRec.cpp
   AMReX_BCRec.H
   AMReX_PhysBCFunct.cpp
   AMReX_PhysBCFunct.H
   AMReX_BCUtil.cpp
   AMReX_BCUtil.H
   AMReX_BC_TYPES.H
   AMReX_FilCC_${DIM}D_C.H
   AMReX_FilCC_C.H
   AMReX_FilCC_C.cpp
   AMReX_FilND_C.H
   AMReX_FilND_C.cpp
   # Plotfile ----------------------------------------------------------------
   AMReX_PlotFileUtil.cpp
   AMReX_PlotFileUtil.H
   AMReX_PlotFileDataImpl.H
   AMReX_PlotFileDataImpl.cpp
   # GPU --------------------------------------------------------------------
   AMReX_Gpu.H
   AMReX_GpuQualifiers.H
   AMReX_GpuPrint.H
   AMReX_GpuAssert.H
   AMReX_GpuTypes.H
   AMReX_GpuControl.H
   AMReX_GpuControl.cpp
   AMReX_GpuFuse.H
   AMReX_GpuFuse.cpp
   AMReX_GpuLaunch.H
   AMReX_GpuLaunch.cpp
   AMReX_GpuLaunchGlobal.H
   AMReX_GpuLaunchMacrosG.H
   AMReX_GpuLaunchMacrosC.H
   AMReX_GpuLaunchFunctsG.H
   AMReX_GpuLaunchFunctsC.H
   AMReX_GpuError.H
   AMReX_GpuDevice.H
   AMReX_GpuDevice.cpp
   AMReX_GpuAtomic.H
   AMReX_GpuUtility.H
   AMReX_GpuUtility.cpp
   AMReX_GpuAsyncArray.H
   AMReX_GpuAsyncArray.cpp
   AMReX_GpuElixir.H
   AMReX_GpuElixir.cpp
   AMReX_GpuMemory.H
   AMReX_GpuRange.H
   AMReX_GpuReduce.H
   AMReX_GpuAllocators.H
   AMReX_GpuContainers.H
   # CUDA --------------------------------------------------------------------
   AMReX_CudaGraph.H
   # Machine model -----------------------------------------------------------
   AMReX_Machine.H
   AMReX_Machine.cpp
   # Memory pool -------------------------------------------------------------
   AMReX_MemPool.cpp
   AMReX_MemPool.H
   # Profiling ---------------------------------------------------------------
   AMReX_BLProfiler.cpp
   AMReX_BLBackTrace.cpp
   )

# Fortran stuff
if (ENABLE_FORTRAN)
   target_sources( amrex
      PRIVATE
      AMReX_filcc_f.H
      AMReX_BLutil_F.F90
      AMReX_BLProfiler_F.F90
      AMReX_FILCC_${DIM}D.F90
      AMReX_filcc_mod.F90
      AMReX_omp_mod.F90
      AMReX_acc_mod.F90
      AMReX_fort_mod.F90
      AMReX_constants_mod.f90
      AMReX_error_mod.F90
      AMReX_parmparse_mod.F90
      AMReX_string_mod.F90
      AMReX_bc_types_mod.F90
      AMReX_ParallelDescriptor_F.F90
      AMReX_io_mod.F90
      AMReX_mempool_mod.F90
      )
endif ()
# LAZY mode not supported in CMake yet
# ifeq ($(LAZY),TRUE)
# C$(AMREX_BASE)_sources += AMReX_Lazy.cpp
# C$(AMREX_BASE)_headers += AMReX_Lazy.H
# endif

# Memory Profiler
if (ENABLE_MEM_PROFILE)
   target_sources(amrex PRIVATE AMReX_MemProfiler.cpp AMReX_MemProfiler.H )
endif ()

# Tiny Profiler
if (ENABLE_TINY_PROFILE)
   target_sources(amrex PRIVATE AMReX_TinyProfiler.cpp AMReX_TinyProfiler.H )
endif ()
