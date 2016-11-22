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

set(CXX_source_files Arena.cpp BArena.cpp BaseFab.cpp BCRec.cpp BLBackTrace.cpp BoxArray.cpp Box.cpp BoxDomain.cpp BoxLib.cpp BoxList.cpp CArena.cpp CoordSys.cpp DistributionMapping.cpp FabArray.cpp FabConv.cpp FArrayBox.cpp FPC.cpp Geometry.cpp MultiFabUtil.cpp IArrayBox.cpp IndexType.cpp IntVect.cpp iMultiFab.cpp MemPool.cpp MultiFab.cpp NFiles.cpp Orientation.cpp ParallelDescriptor.cpp ParmParse.cpp Periodicity.cpp PhysBCFunct.cpp PlotFileUtil.cpp RealBox.cpp UseCount.cpp Utility.cpp VisMF.cpp)

set(F77_source_files BLBoxLib_F.f bl_flush.f BLParmParse_F.f BLutil_F.f)
set(FPP_source_files COORDSYS_${BL_SPACEDIM}D.F FILCC_${BL_SPACEDIM}D.F)
set(F90PP_source_files bl_fort_module.F90)
set(F90_source_files mempool_f.f90 threadbox.f90 MultiFabUtil_${BL_SPACEDIM}d.f90 BaseFab_nd.f90)

set(CXX_header_files Arena.H Array.H ArrayLim.H BArena.H BaseFab.H BCRec.H BC_TYPES.H BLassert.H BLBackTrace.H BLFort.H BLProfiler.H BoxArray.H BoxDomain.H Box.H BoxLib.H BoxList.H CArena.H ccse-mpi.H CONSTANTS.H CoordSys.H DistributionMapping.H FabArray.H FabConv.H FArrayBox.H FPC.H Geometry.H MultiFabUtil.H IArrayBox.H IndexType.H IntVect.H Looping.H iMultiFab.H MemPool.H MultiFab.H NFiles.H Orientation.H ParallelDescriptor.H ParmParse.H PArray.H Periodicity.H PList.H PlotFileUtil.H Pointers.H RealBox.H REAL.H SPACE.H Tuple.H UseCount.H Utility.H VisMF.H winstd.H PhysBCFunct.H)

set(F77_header_files bc_types.fi)
set(FPP_header_files COORDSYS_F.H SPACE_F.H BaseFab_f.H)
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
#   target_link_libraries(boxlib_tBA boxlib ${BOXLIB_EXTRA_LIBRARIES} ${UnitTest_LIBRARIES})
#   add_ccse_test("boxlib_1_tBA" "boxlib_tBA" KIND unit "${tBA_args}")

#   set(tDir_args "Testing/Temporary/dir_created_by_tDir_remove_me")
#   add_executable(boxlib_tDir test/tDir.cpp)
#   target_link_libraries(boxlib_tDir boxlib ${BOXLIB_EXTRA_LIBRARIES} ${UnitTest_LIBRARIES})
#   add_ccse_test("boxlib_2_tDir" "boxlib_tDir" KIND unit "${tDir_args}")

#   set(tDM_args ${BOXLIB_FILE_DIR}/ba.213 ${BOXLIB_FILE_DIR}/ba.60)
#   add_executable(boxlib_tDM test/tDM.cpp)
#   target_link_libraries(boxlib_tDM boxlib ${BOXLIB_EXTRA_LIBRARIES} ${UnitTest_LIBRARIES})
#   add_ccse_test("boxlib_3_tDM" "boxlib_tDM" KIND unit "${tDM_args}")

#   if(BOXLIB_IS_PARALLEL)
#       set(tPD_args)
#       add_executable(boxlib_tPD test/tPD.cpp)
#       target_link_libraries(boxlib_tPD boxlib ${BOXLIB_EXTRA_LIBRARIES} ${UnitTest_LIBRARIES})
#       add_ccse_test("boxlib_4_tPD" "boxlib_tPD" KIND unit NPROCS 8 "${tPD_args}")
#   endif()
#endif()
   
