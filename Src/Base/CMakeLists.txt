# -*- mode: cmake -*-

include(TestManager)

#
# Define a project name
# After this command the following variables are defined
#   CBOXLIB_SOURCE_DIR
#   CBOXLIB_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(CBOXLIB)

include_directories(${CBOXLIB_INCLUDE_DIRS})

set(CXX_source_files AMReX_Arena.cpp AMReX_BArena.cpp AMReX_BaseFab.cpp AMReX_BCRec.cpp AMReX_BLBackTrace.cpp AMReX_BoxArray.cpp AMReX_Box.cpp AMReX_BoxIterator.cpp AMReX_BoxDomain.cpp AMReX.cpp AMReX_BoxList.cpp AMReX_CArena.cpp AMReX_CoordSys.cpp AMReX_DistributionMapping.cpp AMReX_FabArrayBase.cpp AMReX_FabConv.cpp AMReX_FArrayBox.cpp AMReX_FPC.cpp AMReX_Geometry.cpp AMReX_MultiFabUtil.cpp AMReX_IArrayBox.cpp AMReX_IndexType.cpp AMReX_IntVect.cpp AMReX_iMultiFab.cpp AMReX_MemPool.cpp AMReX_MultiFab.cpp AMReX_MFCopyDescriptor.cpp AMReX_MFIter.cpp AMReX_NFiles.cpp AMReX_Orientation.cpp AMReX_ParallelDescriptor.cpp AMReX_ParmParse.cpp AMReX_Periodicity.cpp AMReX_PhysBCFunct.cpp AMReX_PlotFileUtil.cpp AMReX_RealBox.cpp AMReX_RealVect.cpp AMReX_Utility.cpp AMReX_VisMF.cpp AMReX_parmparse_fi.cpp AMReX_error_fi.cpp AMReX_BLProfiler.cpp AMReX_parstream.cpp AMReX_SPMD.cpp AMReX_TinyProfiler.cpp AMReX_BCUtil.cpp)

set(F77_source_files AMReX_BLProfiler_F.f AMReX_BLBoxLib_F.f AMReX_bl_flush.f AMReX_BLParmParse_F.f AMReX_BLutil_F.f)
set(FPP_source_files AMReX_COORDSYS_${BL_SPACEDIM}D.F AMReX_FILCC_${BL_SPACEDIM}D.F)
set(F90PP_source_files AMReX_filcc_mod.F90 AMReX_omp_mod.F90 AMReX_fort_mod.F90 AMReX_parmparse_mod.F90 AMReX_string_mod.F90 AMReX_error_mod.F90 AMReX_bc_types_mod.F90)
set(F90_source_files AMReX_mempool_f.f90 AMReX_MultiFabUtil_${BL_SPACEDIM}d.f90 AMReX_BaseFab_nd.f90)

set(CXX_header_files AMReX_Arena.H AMReX_Array.H AMReX_ArrayLim.H AMReX_BArena.H AMReX_BaseFab.H AMReX_BCRec.H AMReX_BC_TYPES.H AMReX_BLassert.H AMReX_BLBackTrace.H AMReX_BLFort.H AMReX_BLProfiler.H AMReX_BoxArray.H AMReX_BoxDomain.H AMReX_Box.H AMReX_BoxIterator.H AMReX.H AMReX_BoxList.H AMReX_CArena.H AMReX_ccse-mpi.H AMReX_CONSTANTS.H AMReX_CoordSys.H AMReX_DistributionMapping.H AMReX_FabArray.H AMReX_FabArrayBase.H AMReX_FACopyDescriptor.H AMReX_FabConv.H AMReX_FabFactory.H AMReX_FArrayBox.H AMReX_FPC.H AMReX_Geometry.H AMReX_MultiFabUtil.H AMReX_IArrayBox.H AMReX_IndexType.H AMReX_IntVect.H AMReX_Looping.H AMReX_iMultiFab.H AMReX_MakeType.H AMReX_MemPool.H AMReX_MultiFab.H AMReX_MFCopyDescriptor.H AMReX_MFIter.H AMReX_NFiles.H AMReX_Orientation.H AMReX_ParallelDescriptor.H AMReX_ParmParse.H AMReX_Array.H AMReX_Periodicity.H AMReX_PlotFileUtil.H AMReX_Print.H AMReX_RealBox.H AMReX_RealVect.H AMReX_REAL.H AMReX_SPACE.H AMReX_Tuple.H AMReX_TypeTraits.H AMReX_Utility.H AMReX_parstream.H AMReX_VisMF.H AMReX_PhysBCFunct.H AMReX_SPMD.H AMReX_SPMDI.H AMReX_TinyProfiler.H AMReX_BCUtil.H)

set(F77_header_files AMReX_bc_types.fi)
set(FPP_header_files AMReX_COORDSYS_F.H AMReX_SPACE_F.H AMReX_BaseFab_f.H AMReX_filcc_f.H)
set(F90_header_files)

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})
preprocess_boxlib_fortran90(F90PP_out_files ${F90PP_source_files})
set(local_source_files ${CXX_source_files} ${F77_source_files} ${FPP_out_files} ${F90PP_out_files} ${F90_source_files})
set(local_header_files ${CXX_header_files} ${F77_header_files} ${FPP_header_files} ${F90_header_files})
add_library(box_c OBJECT ${local_source_files})

add_install_include_file(${local_header_files})

#if (BUILD_TESTS AND BL_SPACEDIM EQUAL 3)
#   # Add UnitTest include directories
#   include_directories(${UnitTest_INCLUDE_DIRS})

#   set(BOXLIB_FILE_DIR "${BOXLIB_SOURCE_DIR}/test")

#   set(tBA_args ${BOXLIB_FILE_DIR}/ba.213 ${BOXLIB_FILE_DIR}/ba.60)
#   add_executable(boxlib_tBA test/tBA.cpp)
#   target_link_libraries(boxlib_tBA boxlib ${AMREX_EXTRA_LIBRARIES} ${UnitTest_LIBRARIES})
#   add_ccse_test("boxlib_1_tBA" "boxlib_tBA" KIND unit "${tBA_args}")

#   set(tDir_args "Testing/Temporary/dir_created_by_tDir_remove_me")
#   add_executable(boxlib_tDir test/tDir.cpp)
#   target_link_libraries(boxlib_tDir boxlib ${AMREX_EXTRA_LIBRARIES} ${UnitTest_LIBRARIES})
#   add_ccse_test("boxlib_2_tDir" "boxlib_tDir" KIND unit "${tDir_args}")

#   set(tDM_args ${BOXLIB_FILE_DIR}/ba.213 ${BOXLIB_FILE_DIR}/ba.60)
#   add_executable(boxlib_tDM test/tDM.cpp)
#   target_link_libraries(boxlib_tDM boxlib ${AMREX_EXTRA_LIBRARIES} ${UnitTest_LIBRARIES})
#   add_ccse_test("boxlib_3_tDM" "boxlib_tDM" KIND unit "${tDM_args}")

#   if(BOXLIB_IS_PARALLEL)
#       set(tPD_args)
#       add_executable(boxlib_tPD test/tPD.cpp)
#       target_link_libraries(boxlib_tPD boxlib ${AMREX_EXTRA_LIBRARIES} ${UnitTest_LIBRARIES})
#       add_ccse_test("boxlib_4_tPD" "boxlib_tPD" KIND unit NPROCS 8 "${tPD_args}")
#   endif()
#endif()
   
