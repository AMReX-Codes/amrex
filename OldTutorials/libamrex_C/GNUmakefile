AMREX_INSTALL_DIR ?= ../../tmp_install_dir

#CXX = CC
#FC = ftn

CXX = mpicxx
FC = gfortran

CPPFLAGS = -I$(AMREX_INSTALL_DIR)/include
CXXFLAGS = -O2 -std=c++11
FFLAGS = -O2
LDFLAGS = -L$(AMREX_INSTALL_DIR)/lib

LIBRARIES = -lamrex

LIBRARIES += -lgfortran
#LIBRARIES += -lifcore

default: main.exe

main.exe: main.o MyAmr.o myf.o
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS) $(LIBRARIES)

main.o: main.cpp MyAmr.H
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

MyAmr.o: MyAmr.cpp MyAmr.H
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

myf.o: myf.f90
	$(FC) -o $@ -c $(FFLAGS) $(CPPFLAGS) $<

.PHONY: clean realclean

clean:
	$(RM) *.o

realclean: clean
	$(RM) main.exe
