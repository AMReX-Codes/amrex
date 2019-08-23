AMREX_INSTALL_DIR ?= ../../tmp_install_dir

#CXX = CC
#FC = ftn

CXX = mpicxx
FC = gfortran

CPPFLAGS = -I. -I$(AMREX_INSTALL_DIR)/include
CXXFLAGS = -O2 -std=c++11
FFLAGS = -O2
LDFLAGS = -L$(AMREX_INSTALL_DIR)/lib

LIBRARIES = -lamrex

LIBRARIES += -lgfortran
#LIBRARIES += -lifcore

default: main.exe

main.exe: main.o test_parameters.o my_func.o
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS) $(LIBRARIES)

main.o: main.cpp MyParams.H
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

test_parameters.o: test_parameters.cpp MyParams.H
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

my_func.o: my_func.f90
	$(FC) -o $@ -c $(FFLAGS) $(CPPFLAGS) $<

.PHONY: clean realclean

clean:
	$(RM) *.o

realclean: clean
	$(RM) main.exe
