# Install script for directory: /home/runner/work/amrex/amrex

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/runner/work/amrex/amrex/installdir")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/runner/work/amrex/amrex/build_test/Src/cmake_install.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/AMReX" TYPE FILE FILES
    "/home/runner/work/amrex/amrex/build_test/lib/cmake/AMReX/AMReXConfig.cmake"
    "/home/runner/work/amrex/amrex/build_test/lib/cmake/AMReX/AMReXConfigVersion.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/runner/work/amrex/amrex/build_test/Src/libamrex_3d.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_ccse-mpi.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Math.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Algorithm.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Any.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Array.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_BlockMutex.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Enum.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuComplex.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Order.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_SmallMatrix.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_ConstexprFor.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Vector.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_TableData.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Tuple.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_TypeList.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Demangle.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Exception.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Extension.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_PODVector.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_ParmParse.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Functional.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Stack.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_String.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Utility.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_FileSystem.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_ValLocPair.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Reduce.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Scan.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Partition.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Morton.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Random.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_RandomEngine.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_BLassert.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_ArrayLim.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_REAL.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_INT.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_CONSTANTS.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_SPACE.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_DistributionMapping.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_ParallelDescriptor.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_OpenMP.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_ParallelReduce.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_ForkJoin.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_ParallelContext.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_VisMFBuffer.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_VisMF.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_AsyncOut.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_BackgroundThread.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Arena.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_BArena.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_CArena.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_PArena.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_SArena.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_DataAllocator.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_BLProfiler.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_BLBackTrace.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_BLFort.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_NFiles.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_parstream.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_ANSIEscCode.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_FabConv.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_FPC.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_VectorIO.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Print.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_IntConv.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_IOFormat.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Box.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_BoxIterator.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Dim3.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_IntVect.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_IndexType.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Loop.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Loop.nolint.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Orientation.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Periodicity.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_RealBox.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_RealVect.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_BoxList.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_BoxArray.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_BoxDomain.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_FArrayBox.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_IArrayBox.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_BaseFab.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Array4.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_MakeType.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_TypeTraits.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_FabDataType.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_FabFactory.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_BaseFabUtility.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_MultiFab.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_MFCopyDescriptor.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_iMultiFab.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_FabArrayBase.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_MFIter.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_FabArray.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_FACopyDescriptor.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_FabArrayCommI.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_FBI.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_PCI.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_FabArrayUtility.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_LayoutData.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_CoordSys.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_COORDSYS_3D_C.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_COORDSYS_C.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Geometry.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_MultiFabUtil.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_MultiFabUtilI.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_MultiFabUtil_3D_C.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_MultiFabUtil_nd_C.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_MultiFabUtil_C.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_BCRec.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_PhysBCFunct.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_BCUtil.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_BC_TYPES.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_FilCC_3D_C.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_FilCC_C.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_FilFC_3D_C.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_FilFC_C.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_FilND_C.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_NonLocalBC.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_NonLocalBCImpl.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_PlotFileUtil.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_PlotFileDataImpl.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_FEIntegrator.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_IntegratorBase.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_RKIntegrator.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_TimeIntegrator.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_RungeKutta.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Gpu.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuQualifiers.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuKernelInfo.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuPrint.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuAssert.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuTypes.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuControl.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuLaunch.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuLaunch.nolint.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuLaunchGlobal.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuLaunchMacrosG.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuLaunchMacrosG.nolint.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuLaunchMacrosC.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuLaunchMacrosC.nolint.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuLaunchFunctsG.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuLaunchFunctsC.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuLaunchFunctsSIMD.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuError.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuDevice.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuBuffer.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuAtomic.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuUtility.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuAsyncArray.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuElixir.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuMemory.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuRange.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuReduce.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuAllocators.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_GpuContainers.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_MFParallelFor.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_MFParallelForC.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_MFParallelForG.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_SIMD.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_TagParallelFor.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_CTOParallelForImpl.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_ParReduce.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_CudaGraph.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Machine.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_MemPool.H"
    "/home/runner/work/amrex/amrex/Src/Base/Parser/AMReX_Parser.H"
    "/home/runner/work/amrex/amrex/Src/Base/Parser/AMReX_Parser_Exe.H"
    "/home/runner/work/amrex/amrex/Src/Base/Parser/AMReX_Parser_Y.H"
    "/home/runner/work/amrex/amrex/Src/Base/Parser/amrex_parser.lex.nolint.H"
    "/home/runner/work/amrex/amrex/Src/Base/Parser/amrex_parser.tab.nolint.H"
    "/home/runner/work/amrex/amrex/Src/Base/Parser/AMReX_IParser.H"
    "/home/runner/work/amrex/amrex/Src/Base/Parser/AMReX_IParser_Exe.H"
    "/home/runner/work/amrex/amrex/Src/Base/Parser/AMReX_IParser_Y.H"
    "/home/runner/work/amrex/amrex/Src/Base/Parser/amrex_iparser.lex.nolint.H"
    "/home/runner/work/amrex/amrex/Src/Base/Parser/amrex_iparser.tab.nolint.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_LUSolver.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_Slopes_K.H"
    "/home/runner/work/amrex/amrex/Src/Base/AMReX_BaseFwd.H"
    "/home/runner/work/amrex/amrex/Src/Boundary/AMReX_FabSet.H"
    "/home/runner/work/amrex/amrex/Src/Boundary/AMReX_BndryRegister.H"
    "/home/runner/work/amrex/amrex/Src/Boundary/AMReX_Mask.H"
    "/home/runner/work/amrex/amrex/Src/Boundary/AMReX_MultiMask.H"
    "/home/runner/work/amrex/amrex/Src/Boundary/AMReX_BndryData.H"
    "/home/runner/work/amrex/amrex/Src/Boundary/AMReX_BoundCond.H"
    "/home/runner/work/amrex/amrex/Src/Boundary/AMReX_InterpBndryData.H"
    "/home/runner/work/amrex/amrex/Src/Boundary/AMReX_LO_BCTYPES.H"
    "/home/runner/work/amrex/amrex/Src/Boundary/AMReX_InterpBndryData_K.H"
    "/home/runner/work/amrex/amrex/Src/Boundary/AMReX_InterpBndryData_3D_K.H"
    "/home/runner/work/amrex/amrex/Src/Boundary/AMReX_LOUtil_K.H"
    "/home/runner/work/amrex/amrex/Src/Boundary/AMReX_YAFluxRegister.H"
    "/home/runner/work/amrex/amrex/Src/Boundary/AMReX_YAFluxRegister_K.H"
    "/home/runner/work/amrex/amrex/Src/Boundary/AMReX_YAFluxRegister_3D_K.H"
    "/home/runner/work/amrex/amrex/Src/Boundary/AMReX_BoundaryFwd.H"
    "/home/runner/work/amrex/amrex/Src/Boundary/AMReX_EdgeFluxRegister.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_AmrCore.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_Cluster.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_ErrorList.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_FillPatchUtil.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_FillPatchUtil_I.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_FillPatcher.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_FluxRegister.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_InterpBase.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_MFInterpolater.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_Interpolater.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_TagBox.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_AmrMesh.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_FluxReg_3D_C.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_FluxReg_C.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_Interp_C.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_Interp_3D_C.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_MFInterp_C.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_MFInterp_3D_C.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_InterpFaceRegister.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_InterpFaceReg_C.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_InterpFaceReg_3D_C.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_AmrCoreFwd.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_AmrParGDB.H"
    "/home/runner/work/amrex/amrex/Src/AmrCore/AMReX_AmrParticles.H"
    "/home/runner/work/amrex/amrex/Src/Amr/AMReX_LevelBld.H"
    "/home/runner/work/amrex/amrex/Src/Amr/AMReX_Amr.H"
    "/home/runner/work/amrex/amrex/Src/Amr/AMReX_AmrLevel.H"
    "/home/runner/work/amrex/amrex/Src/Amr/AMReX_Derive.H"
    "/home/runner/work/amrex/amrex/Src/Amr/AMReX_StateData.H"
    "/home/runner/work/amrex/amrex/Src/Amr/AMReX_PROB_AMR_F.H"
    "/home/runner/work/amrex/amrex/Src/Amr/AMReX_StateDescriptor.H"
    "/home/runner/work/amrex/amrex/Src/Amr/AMReX_AuxBoundaryData.H"
    "/home/runner/work/amrex/amrex/Src/Amr/AMReX_Extrapolater.H"
    "/home/runner/work/amrex/amrex/Src/Amr/AMReX_extrapolater_K.H"
    "/home/runner/work/amrex/amrex/Src/Amr/AMReX_extrapolater_3D_K.H"
    "/home/runner/work/amrex/amrex/Src/Amr/AMReX_AmrFwd.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLMG.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLMG_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLMG_3D_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLMGBndry.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLLinOp.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLLinOp_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLCellLinOp.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLNodeLinOp.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLNodeLinOp_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLNodeLinOp_3D_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLCellABecLap.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLCellABecLap_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLCellABecLap_3D_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLCGSolver.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_PCGSolver.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLABecLaplacian.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLABecLap_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLABecLap_3D_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLALaplacian.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLALap_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLALap_3D_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLPoisson.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLPoisson_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLPoisson_3D_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/AMReX_GMRES.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/AMReX_GMRES_MLMG.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/AMReX_GMRES_MV.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/AMReX_Smoother_MV.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/AMReX_Algebra.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/AMReX_AlgPartition.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/AMReX_AlgVector.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/AMReX_AlgVecUtil.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/AMReX_CSR.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/AMReX_SpMatrix.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/AMReX_SpMatUtil.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/AMReX_SpMV.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLMG_2D_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLPoisson_2D_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLALap_2D_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLCurlCurl.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLCurlCurl_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLEBNodeFDLaplacian.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLEBNodeFDLap_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLEBNodeFDLap_3D_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLNodeTensorLaplacian.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLNodeTensorLap_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLNodeTensorLap_3D_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLNodeABecLaplacian.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLNodeABecLap_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLNodeABecLap_3D_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLNodeLaplacian.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLNodeLap_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLNodeLap_3D_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLTensorOp.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLTensor_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/MLMG/AMReX_MLTensor_3D_K.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/OpenBC/AMReX_OpenBC.H"
    "/home/runner/work/amrex/amrex/Src/LinearSolvers/OpenBC/AMReX_OpenBC_K.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_Particles.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_ParticleContainer.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_SparseBins.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_ParGDB.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_Particle_mod_K.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_TracerParticles.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_NeighborParticles.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_NeighborParticlesI.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_NeighborList.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_Particle.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_ParticleInit.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_ParticleContainerI.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_ParIter.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_ParticleMPIUtil.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_ParticleUtil.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_StructOfArrays.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_ArrayOfStructs.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_ParticleTile.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_MakeParticle.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_NeighborParticlesCPUImpl.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_NeighborParticlesGPUImpl.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_ParticleBufferMap.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_ParticleCommunication.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_ParticleInterpolators.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_ParticleReduce.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_ParticleMesh.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_ParticleLocator.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_ParticleIO.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_DenseBins.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_BinIterator.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_ParticleTransformation.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_WriteBinaryParticleData.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_ParticleContainerBase.H"
    "/home/runner/work/amrex/amrex/Src/Particle/AMReX_ParticleArray.H"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/AMReX/AMReXTargets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/AMReX/AMReXTargets.cmake"
         "/home/runner/work/amrex/amrex/build_test/CMakeFiles/Export/2260e541ece776bcef17e59de6c71ec8/AMReXTargets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/AMReX/AMReXTargets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/AMReX/AMReXTargets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/AMReX" TYPE FILE FILES "/home/runner/work/amrex/amrex/build_test/CMakeFiles/Export/2260e541ece776bcef17e59de6c71ec8/AMReXTargets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/AMReX" TYPE FILE FILES "/home/runner/work/amrex/amrex/build_test/CMakeFiles/Export/2260e541ece776bcef17e59de6c71ec8/AMReXTargets-release.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(CREATE_LINK
           libamrex_3d.a
           "/home/runner/work/amrex/amrex/installdir/lib/libamrex.a"
           COPY_ON_ERROR SYMBOLIC)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/amrex" TYPE DIRECTORY FILES
    "/home/runner/work/amrex/amrex/Tools/C_scripts"
    "/home/runner/work/amrex/amrex/Tools/typechecker"
    USE_SOURCE_PERMISSIONS)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/AMReX/AMReXCMakeModules" TYPE DIRECTORY FILES "/home/runner/work/amrex/amrex/Tools/CMake/" USE_SOURCE_PERMISSIONS)
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/runner/work/amrex/amrex/build_test/Tests/cmake_install.cmake")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/runner/work/amrex/amrex/build_test/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
if(CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_COMPONENT MATCHES "^[a-zA-Z0-9_.+-]+$")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
  else()
    string(MD5 CMAKE_INST_COMP_HASH "${CMAKE_INSTALL_COMPONENT}")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INST_COMP_HASH}.txt")
    unset(CMAKE_INST_COMP_HASH)
  endif()
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/runner/work/amrex/amrex/build_test/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
