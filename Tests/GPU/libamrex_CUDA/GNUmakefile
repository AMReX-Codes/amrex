AMREX_INSTALL_DIR ?= ../../../tmp_install_dir
CUDA_ARCH ?= 60

CXX = nvcc

includes = -I$(AMREX_INSTALL_DIR)/include
CXXFLAGS = -ccbin=g++ -Xcompiler='-g -O3 -std=c++14' --std=c++14 -Wno-deprecated-gpu-targets -m64 -arch=compute_$(CUDA_ARCH) -code=sm_$(CUDA_ARCH) -maxrregcount=255 -lineinfo --ptxas-options=-O3,-v --use_fast_math -dc -x cu --expt-relaxed-constexpr --expt-extended-lambda
LINKFLAGS = -Wno-deprecated-gpu-targets -m64 -arch=compute_$(CUDA_ARCH) -code=sm_$(CUDA_ARCH) -maxrregcount=255 -lineinfo --ptxas-options=-O3,-v --use_fast_math
LDFLAGS = -L$(AMREX_INSTALL_DIR)/lib

LIBRARIES = -lamrex -lgfortran -lcuda

# libamrex.a may contain unresolved symbols, amrex_probinit and getLevelBld()
# -Xlinker='--unresolved-symbols=ignore-in-object-files'
LINKFLAGS += -Xlinker='--warn-unresolved-symbols'

default: main.exe

main.exe: main.o
	$(CXX) -o $@ $< $(includes) $(LINKFLAGS) $(LDFLAGS) $(LIBRARIES)

main.o: main.cpp
	$(CXX) -o $@ -c $(CXXFLAGS) $(includes) $<

.PHONY: clean realclean

clean:
	$(RM) *.o

realclean: clean
	$(RM) main.exe
