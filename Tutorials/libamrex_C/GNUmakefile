AMREX_INSTALL_DIR ?= ../../tmp_install_dir

default: main.exe

main.exe: main.o MyAmr.o myf.o
	mpicxx -o $@ $^ -L$(AMREX_INSTALL_DIR)/lib -lamrex -lgfortran

main.o: main.cpp MyAmr.H
	mpicxx -o $@ -c -O2 -std=c++11 -I$(AMREX_INSTALL_DIR)/include $<

MyAmr.o: MyAmr.cpp MyAmr.H
	mpicxx -o $@ -c -O2 -std=c++11 -I$(AMREX_INSTALL_DIR)/include $<

myf.o: myf.f90
	gfortran -o $@ -c -O2 -I$(AMREX_INSTALL_DIR)/include $<

.PHONY: clean realclean

clean:
	$(RM) *.o

realclean: clean
	$(RM) main.exe
