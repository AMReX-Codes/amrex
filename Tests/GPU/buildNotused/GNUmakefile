############################

objects = main.o MyKernel.o MyKernelB.o Notused.o

bsrcs = Notused.o
srcs = MyKernel.o MyKernelB.o 
main = main.o

app = link.app
lib = lib.a

#Defaults, used to make with CUDA as an executable.

CXX = nvcc
CXXFLAGS = -ccbin=g++ -Xcompiler='-g -O3 -std=c++11 --std=c++11' --std=c++11 -m64 -arch=compute_70 -code=sm_70 -dc -x cu -DUSE_CUDA 

USE_CUDA = -DUSE_CUDA

############################
# SIMPLE, MAKE AN EXEC TARGETS
# 
# all USES nvcc w/ CUDA
# gcc USES g++  w/o CUDA
#
# broke & gccbrokei
# show how compiling fails when
#    making an executable including
#    the object with an undefined file
# ########################### 
 
all: $(objects)
	nvcc -m64 -arch=compute_70 -code=sm_70 -lcuda -lnvidia-ml $(srcs) $(main) -o $(app) 

broke: $(objects)
	nvcc -m64 -arch=compute_70 -code=sm_70 -lcuda -lnvidia-ml $(objects) -o $(app) 

# Make w/o CUDA
gcc: CXX=g++
gcc: CXXFLAGS= -c
gcc: $(objects)
	g++ $(srcs) $(main) -o $(app) 

gccbroke: CXX=g++
gccbroke: CXXFLAGS= -c
gccbroke: $(objects)
	g++ $(objects) -o $(app) 

# ============================
# Building without CUDA
# Compiles correctly w/ g++
# Note: Notused included in library

gccexec: CXX=g++
gccexec: CXXFLAGS= -c
gccexec: gcclibrary $(main)
	g++ -lstdc++ -latomic $(main) $(lib) -o $(app).gpp

gcclibrary: CXX=g++
gcclibrary: CXXFLAGS= -c
gcclibrary: $(srcs) $(bsrcs)
	$(SILENT) $(AR) qc $(lib) $^
	$(SILENT) ranlib $(lib)

# ============================
# Building without CUDA
# Compiles correctly w/ nvcc
# Note: Notused included in library

nvccexec: COMP = nvcc
nvccexec: CXXFLAGS= -ccbin=g++ -Xcompiler='-g -O3 -std=c++11 --std=c++11' --std=c++11 -m64 -c
nvccexec: nvcclibrary $(main)
	nvcc -lstdc++ -latomic $(main) $(lib) -o $(app).gpp

nvcclibrary: CXXFLAGS= -ccbin=g++ -Xcompiler='-g -O3 -std=c++11 --std=c++11' --std=c++11 -m64 -c
nvcclibrary: $(srcs) $(bsrcs)
	$(SILENT) $(AR) qc $(lib) $^
	$(SILENT) ranlib $(lib)

# g++ build with CUDA
# requires dlink of library
#   and main.o, which forces
#   definition of notdefined
#   Leads to failure
# =========================

cudaexec: cudalibrary $(main)
	nvcc -ccbin=g++ -Xcompiler=-fPIC --std=c++11 -m64 -arch=compute_70 -code=sm_70 -dlink main.o $(lib) -o link.o
	g++ -gopt -m64 -Wl,--unresolved-symbols=ignore-in-object-files -I. -I/sw/summit/cuda/9.2.148/include -L. -L/sw/summit/cuda/9.2.148/lib64 -L/usr/lib64/nvidia -lcuda -lcudart -lcudadevrt -lnvidia-ml -lstdc++ -latomic $(main) $(bsrcs) $(lib) link.o -o $(app).gpp

cudaexec: CXX = nvcc
cudaexec: CXXFLAGS = -ccbin=g++ -Xcompiler='-g -O3 -std=c++11 --std=c++11' --std=c++11 -m64 -arch=compute_70 -code=sm_70 -dc -x cu -DUSE_CUDA 
cudalibrary: $(bsrcs) $(srcs)
	nvcc -ccbin=g++ -Xcompiler='-g -O3 -std=c++11 --std=c++11' --std=c++11 -m64 -arch=compute_70 -code=sm_70 --lib -lcuda -lnvidia-ml $(bsrcs) $(srcs) -o $(lib)

# =========================

# MAKE INDIVIDUAL OBJECT FILES
# ==============================

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ 

clean:
	rm -f *.o *.mod *.a $(app)*

realclean: clean
