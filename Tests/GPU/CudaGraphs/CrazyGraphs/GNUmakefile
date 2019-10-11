MAKE_FLAGS := -I$(CUDA_DIR)/include -L$(CUDA_DIR)/lib64

all: main.cpp
	nvcc main.cpp -ccbin=g++ -Xcompiler=' -g -O3  --std=c++14' --std=c++14 --expt-extended-lambda -Wno-deprecated-gpu-targets -m64 -arch=compute_70 -code=sm_70 -x cu -I. $(MAKE_FLAGS) -L. -L/opt/gcc/7.3.0/snos/lib/gcc/x86_64-suse-linux/7.3.0/../../../../lib64/ -lcuda -o graphs.ex 

clean:
	rm graphs.ex

realclean:
	rm graphs.ex
