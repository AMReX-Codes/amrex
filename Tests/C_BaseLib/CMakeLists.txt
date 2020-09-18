
#
# Disabled for NOW
#
return()


if (NOT CMAKE_Fortran_COMPILER_LOADED)
   return()
endif ()

set(_sources AMRProfTestBL.cpp tBA.cpp tDM.cpp tFillFab.cpp tParmParse.cpp  TPROFILER_F.H  tread.cpp  tVisMF.cpp )
list(APPEND _sources AMRPROFTEST_F.H  tCArena.cpp  tFAC.cpp tMFcopy.cpp  tProfiler.cpp tRABcast.cpp tUMap.cpp    )
list(APPEND _sources fillfab.f  t8BIT.cpp tDir.cpp tFB.cpp  tMF.cpp TPROFILER.F tRan.cpp tVisMF2.cpp )

set(_input_files )

setup_test(_sources _input_files)

unset(_sources)
unset(_input_files)
