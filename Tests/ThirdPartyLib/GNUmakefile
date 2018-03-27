# Here we assume gcc and mpich.

# We assume that libamrex.a has been built and installed at
AMREX_DIR := ../../tmp_install_dir

default: main.ex

main.ex: main.o libfoo_cxx.a libbar_f.a
	mpicc -o $@ -O $^ -L. -L$(AMREX_DIR)/lib -lfoo_cxx -lbar_f -lamrex -lmpichf90 -lgfortran -lstdc++ -lm

main.o: main.c
	mpicc -o $@ -c -O $^

libfoo_cxx.a: foo.o
	ar rv $@ $^
	ranlib $@

foo.o: foo.cpp
	mpicxx -o $@ -c -O -std=c++11 $^ -I$(AMREX_DIR)/include

libbar_f.a: bar.o
	ar rv $@ $^
	ranlib $@

bar.o: bar.F90
	mpif90 -o $@ -c -O $^ -I$(AMREX_DIR)/include

clean::
	$(RM) main.ex *.o *.a

.PHONY: clean
