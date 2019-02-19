dimension = get_option('dimension')
if dimension == '3'
    src_dim = [
        'Amr/AMReX_extrapolater_3d.f90',
        'AmrCore/AMReX_FillPatchUtil_3d.F90',
        'AmrCore/AMReX_INTERP_3D.F90',
        'Base/AMReX_FILCC_3D.F90',
    ]
elif dimension == '2'
    src_dim = [
        'Amr/AMReX_extrapolater_2d.f90',
        'AmrCore/AMReX_FillPatchUtil_2d.F90',
        'AmrCore/AMReX_INTERP_2D.F90',
        'Base/AMReX_FILCC_2D.F90',
    ]
elif dimension == '1'
    src_dim = [
        'Amr/AMReX_extrapolater_1d.f90',
        'AmrCore/AMReX_FillPatchUtil_1d.F90',
        'AmrCore/AMReX_INTERP_1D.F90',
        'Base/AMReX_FILCC_1D.F90',
    ]
else
    error('Dimension must be one of [1, 2, 3]; invalid value: '+ get_option('dimension'))
endif

srcs = src_dim + [
    'Amr/AMReX_Amr.cpp',
    'Amr/AMReX_AmrLevel.cpp',
    'Amr/AMReX_AuxBoundaryData.cpp',
    'Amr/AMReX_Derive.cpp',
    'Amr/AMReX_Extrapolater.cpp',
    'Amr/AMReX_StateData.cpp',
    'Amr/AMReX_StateDescriptor.cpp',
    'AmrCore/AMReX_AmrCore.cpp',
    'AmrCore/AMReX_AmrMesh.cpp',
    'AmrCore/AMReX_Cluster.cpp',
    'AmrCore/AMReX_ErrorList.cpp',
    'AmrCore/AMReX_FLUXREG_nd.F90',
    'AmrCore/AMReX_FillPatchUtil.cpp',
    'AmrCore/AMReX_FluxRegister.cpp',
    'AmrCore/AMReX_Interpolater.cpp',
    'AmrCore/AMReX_TagBox.cpp',
    'Base/AMReX.cpp',
    'Base/AMReX_Arena.cpp',
    'Base/AMReX_BArena.cpp',
    'Base/AMReX_BCRec.cpp',
    'Base/AMReX_BCUtil.cpp',
    'Base/AMReX_BLBackTrace.cpp',
    'Base/AMReX_BLBoxLib_F.F90',
    'Base/AMReX_BLProfiler.cpp',
    'Base/AMReX_BLProfiler_F.F90',
    'Base/AMReX_BLutil_F.F90',
    'Base/AMReX_BaseFab.cpp',
    'Base/AMReX_Box.cpp',
    'Base/AMReX_BoxArray.cpp',
    'Base/AMReX_BoxDomain.cpp',
    'Base/AMReX_BoxIterator.cpp',
    'Base/AMReX_BoxList.cpp',
    'Base/AMReX_CArena.cpp',
    'Base/AMReX_CoordSys.cpp',
    'Base/AMReX_CudaAllocators.cpp',
    'Base/AMReX_CudaAsyncArray.cpp',
    'Base/AMReX_CudaAsyncFab.cpp',
    'Base/AMReX_CudaAsyncFabImpl.cpp',
    'Base/AMReX_CudaDevice.cpp',
    'Base/AMReX_CudaElixir.cpp',
    'Base/AMReX_CudaLaunch.cpp',
    'Base/AMReX_CudaUtility.cpp',
    'Base/AMReX_DistributionMapping.cpp',
    'Base/AMReX_FArrayBox.cpp',
    'Base/AMReX_FPC.cpp',
    'Base/AMReX_FabArrayBase.cpp',
    'Base/AMReX_FabConv.cpp',
    'Base/AMReX_ForkJoin.cpp',
    'Base/AMReX_Geometry.cpp',
    'Base/AMReX_GpuControl.cpp',
    'Base/AMReX_IArrayBox.cpp',
    'Base/AMReX_IndexType.cpp',
    'Base/AMReX_IntConv.cpp',
    'Base/AMReX_IntVect.cpp',
    'Base/AMReX_MFCopyDescriptor.cpp',
    'Base/AMReX_MFIter.cpp',
    'Base/AMReX_Machine.cpp',
    'Base/AMReX_MemPool.cpp',
    'Base/AMReX_MultiFab.cpp',
    'Base/AMReX_MultiFabUtil.cpp',
    'Base/AMReX_NFiles.cpp',
    'Base/AMReX_Orientation.cpp',
    'Base/AMReX_ParallelContext.cpp',
    'Base/AMReX_ParallelDescriptor.cpp',
    'Base/AMReX_ParallelDescriptor_F.F90',
    'Base/AMReX_ParmParse.cpp',
    'Base/AMReX_Periodicity.cpp',
    'Base/AMReX_PhysBCFunct.cpp',
    'Base/AMReX_PlotFileDataImpl.cpp',
    'Base/AMReX_PlotFileUtil.cpp',
    'Base/AMReX_RealBox.cpp',
    'Base/AMReX_RealVect.cpp',
    'Base/AMReX_Utility.cpp',
    'Base/AMReX_VectorIO.cpp',
    'Base/AMReX_VisMF.cpp',
    'Base/AMReX_acc_mod.F90',
    'Base/AMReX_bc_types_mod.F90',
    'Base/AMReX_constants_mod.f90',
    'Base/AMReX_error_fi.cpp',
    'Base/AMReX_error_mod.F90',
    'Base/AMReX_filcc_mod.F90',
    'Base/AMReX_fort_mod.F90',
    'Base/AMReX_iMultiFab.cpp',
    'Base/AMReX_io_mod.F90',
    'Base/AMReX_mempool_mod.F90',
    'Base/AMReX_omp_mod.F90',
    'Base/AMReX_parmparse_fi.cpp',
    'Base/AMReX_parmparse_mod.F90',
    'Base/AMReX_parstream.cpp',
    'Base/AMReX_string_mod.F90',
    'Boundary/AMReX_BndryData.cpp',
    'Boundary/AMReX_BndryRegister.cpp',
    'Boundary/AMReX_FabSet.cpp',
    'Boundary/AMReX_InterpBndryData.cpp',
    'Boundary/AMReX_LO_UTIL.F90',
    'Boundary/AMReX_MacBndry.cpp',
    'Boundary/AMReX_Mask.cpp',
    'Boundary/AMReX_MultiMask.cpp',
    'Boundary/AMReX_YAFluxRegister.cpp',
    'Boundary/AMReX_lo_bctypes_mod.F90',
]

mpi = get_option('mpi')
mpi_deps = []
if mpi
    mpi_deps += [
        dependency('mpi', language: 'cpp'),
        dependency('mpi', language: 'fortran'),
    ]
endif

incdirs_amrex = include_directories(
    'Amr',
    'AmrCore',
    'AmrTask/Amr',
    'AmrTask/graph',
    'AmrTask/rts_impls/Serial',
    'Base',
    'Boundary',
    'F_BaseLib',
)

fortran_args = [
    '-DAMREX_SPACEDIM=' + dimension,
    '-cpp',
    '-fno-range-check',
]
if mpi
    fortran_args += [
        '-DBL_USE_MPI',
        '-DAMREX_USE_MPI',
    ]
endif

cpp_args = [
    '-DAMREX_SPACEDIM=' + dimension,
    '-DBL_FORT_USE_LOWERCASE=1',
]
base = static_library(
    'amrex',
    srcs,
    include_directories: incdirs_amrex,
    fortran_args: fortran_args,
    cpp_args: cpp_args,
    dependencies: mpi_deps,
    install: true,
)
