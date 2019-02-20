project('amrex', 'cpp', 'fortran', 'c',
        default_options: ['c_std=c11', 'cpp_std=c++11'])

dimension = get_option('dimension')
if dimension == '3'
    src_dim = [
        'Src/Amr/AMReX_extrapolater_3d.f90',
        'Src/AmrCore/AMReX_FillPatchUtil_3d.F90',
        'Src/AmrCore/AMReX_INTERP_3D.F90',
        'Src/Base/AMReX_FILCC_3D.F90',
    ]
elif dimension == '2'
    src_dim = [
        'Src/Amr/AMReX_extrapolater_2d.f90',
        'Src/AmrCore/AMReX_FillPatchUtil_2d.F90',
        'Src/AmrCore/AMReX_INTERP_2D.F90',
        'Src/Base/AMReX_FILCC_2D.F90',
    ]
elif dimension == '1'
    src_dim = [
        'Src/Amr/AMReX_extrapolater_1d.f90',
        'Src/AmrCore/AMReX_FillPatchUtil_1d.F90',
        'Src/AmrCore/AMReX_INTERP_1D.F90',
        'Src/Base/AMReX_FILCC_1D.F90',
    ]
else
    error('Dimension must be one of [1, 2, 3]; invalid value: '+ get_option('dimension'))
endif

srcs = src_dim + [
    'Src/Amr/AMReX_Amr.cpp',
    'Src/Amr/AMReX_AmrLevel.cpp',
    'Src/Amr/AMReX_AuxBoundaryData.cpp',
    'Src/Amr/AMReX_Derive.cpp',
    'Src/Amr/AMReX_Extrapolater.cpp',
    'Src/Amr/AMReX_StateData.cpp',
    'Src/Amr/AMReX_StateDescriptor.cpp',
    'Src/AmrCore/AMReX_AmrCore.cpp',
    'Src/AmrCore/AMReX_AmrMesh.cpp',
    'Src/AmrCore/AMReX_Cluster.cpp',
    'Src/AmrCore/AMReX_ErrorList.cpp',
    'Src/AmrCore/AMReX_FLUXREG_nd.F90',
    'Src/AmrCore/AMReX_FillPatchUtil.cpp',
    'Src/AmrCore/AMReX_FluxRegister.cpp',
    'Src/AmrCore/AMReX_Interpolater.cpp',
    'Src/AmrCore/AMReX_TagBox.cpp',
    'Src/Base/AMReX.cpp',
    'Src/Base/AMReX_Arena.cpp',
    'Src/Base/AMReX_BArena.cpp',
    'Src/Base/AMReX_BCRec.cpp',
    'Src/Base/AMReX_BCUtil.cpp',
    'Src/Base/AMReX_BLBackTrace.cpp',
    'Src/Base/AMReX_BLBoxLib_F.F90',
    'Src/Base/AMReX_BLProfiler.cpp',
    'Src/Base/AMReX_BLProfiler_F.F90',
    'Src/Base/AMReX_BLutil_F.F90',
    'Src/Base/AMReX_BaseFab.cpp',
    'Src/Base/AMReX_Box.cpp',
    'Src/Base/AMReX_BoxArray.cpp',
    'Src/Base/AMReX_BoxDomain.cpp',
    'Src/Base/AMReX_BoxIterator.cpp',
    'Src/Base/AMReX_BoxList.cpp',
    'Src/Base/AMReX_CArena.cpp',
    'Src/Base/AMReX_CoordSys.cpp',
    'Src/Base/AMReX_CudaAllocators.cpp',
    'Src/Base/AMReX_CudaAsyncArray.cpp',
    'Src/Base/AMReX_CudaAsyncFab.cpp',
    'Src/Base/AMReX_CudaAsyncFabImpl.cpp',
    'Src/Base/AMReX_CudaDevice.cpp',
    'Src/Base/AMReX_CudaElixir.cpp',
    'Src/Base/AMReX_CudaLaunch.cpp',
    'Src/Base/AMReX_CudaUtility.cpp',
    'Src/Base/AMReX_DistributionMapping.cpp',
    'Src/Base/AMReX_FArrayBox.cpp',
    'Src/Base/AMReX_FPC.cpp',
    'Src/Base/AMReX_FabArrayBase.cpp',
    'Src/Base/AMReX_FabConv.cpp',
    'Src/Base/AMReX_ForkJoin.cpp',
    'Src/Base/AMReX_Geometry.cpp',
    'Src/Base/AMReX_GpuControl.cpp',
    'Src/Base/AMReX_IArrayBox.cpp',
    'Src/Base/AMReX_IndexType.cpp',
    'Src/Base/AMReX_IntConv.cpp',
    'Src/Base/AMReX_IntVect.cpp',
    'Src/Base/AMReX_MFCopyDescriptor.cpp',
    'Src/Base/AMReX_MFIter.cpp',
    'Src/Base/AMReX_Machine.cpp',
    'Src/Base/AMReX_MemPool.cpp',
    'Src/Base/AMReX_MultiFab.cpp',
    'Src/Base/AMReX_MultiFabUtil.cpp',
    'Src/Base/AMReX_NFiles.cpp',
    'Src/Base/AMReX_Orientation.cpp',
    'Src/Base/AMReX_ParallelContext.cpp',
    'Src/Base/AMReX_ParallelDescriptor.cpp',
    'Src/Base/AMReX_ParallelDescriptor_F.F90',
    'Src/Base/AMReX_ParmParse.cpp',
    'Src/Base/AMReX_Periodicity.cpp',
    'Src/Base/AMReX_PhysBCFunct.cpp',
    'Src/Base/AMReX_PlotFileDataImpl.cpp',
    'Src/Base/AMReX_PlotFileUtil.cpp',
    'Src/Base/AMReX_RealBox.cpp',
    'Src/Base/AMReX_RealVect.cpp',
    'Src/Base/AMReX_Utility.cpp',
    'Src/Base/AMReX_VectorIO.cpp',
    'Src/Base/AMReX_VisMF.cpp',
    'Src/Base/AMReX_acc_mod.F90',
    'Src/Base/AMReX_bc_types_mod.F90',
    'Src/Base/AMReX_constants_mod.f90',
    'Src/Base/AMReX_error_fi.cpp',
    'Src/Base/AMReX_error_mod.F90',
    'Src/Base/AMReX_filcc_mod.F90',
    'Src/Base/AMReX_fort_mod.F90',
    'Src/Base/AMReX_iMultiFab.cpp',
    'Src/Base/AMReX_io_mod.F90',
    'Src/Base/AMReX_mempool_mod.F90',
    'Src/Base/AMReX_omp_mod.F90',
    'Src/Base/AMReX_parmparse_fi.cpp',
    'Src/Base/AMReX_parmparse_mod.F90',
    'Src/Base/AMReX_parstream.cpp',
    'Src/Base/AMReX_string_mod.F90',
    'Src/Boundary/AMReX_BndryData.cpp',
    'Src/Boundary/AMReX_BndryRegister.cpp',
    'Src/Boundary/AMReX_FabSet.cpp',
    'Src/Boundary/AMReX_InterpBndryData.cpp',
    'Src/Boundary/AMReX_LO_UTIL.F90',
    'Src/Boundary/AMReX_MacBndry.cpp',
    'Src/Boundary/AMReX_Mask.cpp',
    'Src/Boundary/AMReX_MultiMask.cpp',
    'Src/Boundary/AMReX_YAFluxRegister.cpp',
    'Src/Boundary/AMReX_lo_bctypes_mod.F90',
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
    'Src/Amr',
    'Src/AmrCore',
    'Src/AmrTask/Amr',
    'Src/AmrTask/graph',
    'Src/AmrTask/rts_impls/Serial',
    'Src/Base',
    'Src/Boundary',
    'Src/F_BaseLib',
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

libamrex = static_library(
    'amrex',
    srcs,
    include_directories: incdirs_amrex,
    fortran_args: fortran_args,
    cpp_args: cpp_args,
    dependencies: mpi_deps,
    install: true,
)

libamrex_dep = declare_dependency(
    include_directories : incdirs_amrex,
    link_with : libamrex,
)
