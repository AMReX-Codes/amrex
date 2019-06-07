ifeq ($(CUDA_VERSION),10.0)
	MAKE_FLAGS := -I/global/common/cori/software/cuda/10.0/include -L/global/common/cori/software/cuda/10.0/lib64
else ifeq ($(CUDA_VERSION),10.1)
	MAKE_FLAGS := -I/global/common/cori/software/cuda/10.1/include -L/global/common/cori/software/cuda/10.1/lib64 -DCRAZY
endif

all: main.cpp
	nvcc main.cpp -ccbin=g++ -Xcompiler=' -g -O3  --std=c++14' --std=c++14 --expt-extended-lambda -Wno-deprecated-gpu-targets -m64 -arch=compute_70 -code=sm_70 -x cu -I. $(MAKE_FLAGS) -L. -L/opt/gcc/7.3.0/snos/lib/gcc/x86_64-suse-linux/7.3.0/../../../../lib64/ -lcuda -o graphs.ex 

clean:
	rm graphs.ex

realclean:
	rm graphs.ex
