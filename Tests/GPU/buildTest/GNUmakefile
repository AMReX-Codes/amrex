objects = main.o MyKernel.o MyKernelB.o

srcs = MyKernel.o MyKernelB.o
main = main.o

app = link.app
lib = lib

# SIMPLE, MAKE AN EXEC TARGETS
# all USES nvcc
# pgi USES pgfortran
# ============================
 
all: $(objects)
	nvcc -m64 -arch=compute_70 -code=sm_70 -lcuda -lnvidia-ml $(objects) -o $(app) 

pgi: $(objects) 
	pgfortran -gopt -Mnomain -Mcuda=cc70,charstring -noacc -L/usr/lib64/nvidia -lcuda -lnvidia-ml -lstdc++ -latomic $(objects) -o $(app).pgi

# NVCC LIBRARY TARGETS 
# ====================

nvccexec: nvcclibrary $(main)
	nvcc -m64 -arch=compute_70 -code=sm_70 -L/usr/lib64/nvidia -lcuda -lnvidia-ml $(lib).a $(main) -o $(app) 

nvcclibrary: $(srcs)
	nvcc --lib -m64 -arch=compute_70 -code=sm_70 -L/usr/lib64/nvidia -lcuda -lnvidia-ml $(srcs) -o $(lib).a

# g++ LIBRARY TARGETS
# =========================

gccexec: gcclibrary $(main)
	nvcc -ccbin=g++ -Xcompiler=-fPIC --std=c++11 -m64 -arch=compute_70 -code=sm_70 -dlink main.o $(lib).gcc.a -o link.o
	g++ -gopt -m64 -I. -I/sw/summit/cuda/9.1.85/include -L. -L/sw/summit/cuda/9.1.85/lib64 -L/usr/lib64/nvidia -lcuda -lcudart -lcudadevrt -lnvidia-ml -lstdc++ -latomic $(main) $(lib).gcc.a link.o -o $(app).gpp

gcclibrary: $(srcs)
	nvcc -ccbin=g++ -Xcompiler='-g -O3 -std=c++11 --std=c++11' --std=c++11 -m64 -arch=compute_70 -code=sm_70 --lib -lcuda -lnvidia-ml $(srcs) -o $(lib).gcc.a

# pgfortran LIBRARY TARGETS
# =========================

pgiexec: pgilibrary $(main)
	# Requires additional dlink step to link device code before pgfortran link into executable. 
	nvcc -ccbin=g++ -Xcompiler='-g -O3 -std=c++11 --std=c++11' --std=c++11 -m64 -arch=compute_70 -code=sm_70 -dlink main.o $(lib).pgi.a -o link.o 
	pgfortran -gopt -Mnomain -Mcuda=cc70,charstring -Mdclchk -noacc -L/usr/lib64/nvidia -lcuda -lcudart -lcudadevrt -lnvidia-ml -lstdc++ -latomic $(main) $(lib).pgi.a link.o -o $(app).pgi

pgilibrary: $(srcs)
	nvcc -ccbin=g++ -Xcompiler='-g -O3 -std=c++11 --std=c++11' --lib -m64 -arch=compute_70 -code=sm_70 -lcuda -lnvidia-ml $(srcs) -o $(lib).pgi.a


# MAKE OBJECT FILES USING EITHER -dc OR -dlink
# ============================================

%.o: %.cpp
	nvcc -ccbin=g++ -Xcompiler='-g -O3 -std=c++11 --std=c++11' --std=c++11 -m64 -arch=compute_70 -code=sm_70 -dc -x cu $< -o $@ 

clean:
	rm -f *.o *.a $(lib)* $(app)*
