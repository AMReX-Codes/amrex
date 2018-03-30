cxxsources += AMReX_MemPool.cpp
cxxsources += AMReX_CArena.cpp
cxxsources += AMReX_Arena.cpp

f90sources += AMReX_mempool_f.f90

ifdef CUDA
  F90sources += AMReX_CUDA.F90
endif

f90sources += AMReX_constants_mod.f90
