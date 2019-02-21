objects = main.o MyKernel.o MyKernelB.o MyKernel_F.o

srcs = MyKernel.o MyKernelB.o
fsrcs = MyKernel_F.o
main = main.o

app = link.app
lib = lib.a

# SIMPLE, MAKE AN EXEC TARGET
# USED TO CHECK CODE IS CORRECT
# WHEN USING FORTRAN, CANNOT USE NVCC
# ===================================
 
all: $(objects) 
	pgfortran -gopt -Mnomain -Mcuda=cc70,ptxinfo,fastmath,charstring -Mdclchk -noacc -L/usr/lib64/nvidia -lcuda -lnvidia-ml -lstdc++ -latomic $(objects) -o $(app).pgi

# pgfortran LIBRARY TARGETS
# pgilibrary BUILDS THE LIBRARY 
# pgiexec BUILDS AN EXECUTABLE USING CUDA LIBRARY 
# =================================================

exec: library $(main) $(fsrcs)
	# No -dlink step needed with Fortran. pgfortran performs the link operation. 
	pgfortran -gopt -Mnomain -Mcuda=cc70 -Mdclchk -noacc -L/usr/lib64/nvidia -lcuda -lcudart -lcudadevrt -lnvidia-ml -lstdc++ -latomic $(fsrcs) $(main) $(lib) -o $(app).lib

library: $(srcs)
	$(SILENT) $(AR) qc $(lib) $^
	$(SILENT) ranlib $(lib)
	# Build library using nvcc instead of ar & ranlib
	#nvcc -ccbin=g++ -Xcompiler='-g -O3 -std=c++11 --std=c++11' --lib -m64 -arch=compute_70 -code=sm_70 -I. -I/sw/summit/cuda/9.1.85/include -L. -L/sw/summit/cuda/9.1.85/lib64 -L/usr/lib64/nvidia -lcuda -lnvidia-ml $(srcs) $(fsrcs) -o $(lib).pgi.a


# MAKE INDIVIDUAL OBJECT FILES, CPP & F90 
# =======================================

%.o: %.cpp
	nvcc -ccbin=g++ -Xcompiler='-g -O3 -std=c++11 --std=c++11' --std=c++11 -m64 -arch=compute_70 -code=sm_70 -dc -x cu $< -o $@ 

%.o: %.F90
	pgfortran -r -gopt -Mnomain -Mcuda=cc70 -Mdclchk -noacc -c $< -o $@

clean:
	rm -f *.o *.a *.mod $(app) $(lib)
