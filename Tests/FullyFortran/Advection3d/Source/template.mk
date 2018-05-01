# Template for the GNU Compiler Collection on Linux systems
#
# Typical use with mkmf
# mkmf -t linux-gnu.mk -c"-Duse_libMPI -Duse_netCDF" path_names /usr/local/include
AMREX_DIR := ../../../../tmp_install_dir
############
# Commands Macors
FC = /opt/local/bin/mpif90
CC = /opt/local/bin/mpicc
CXX = /opt/local/bin/mpicxx
LD = /opt/local/bin/mpicc $(MAIN_PROGRAM)
FFLAGS +=-I$(AMREX_DIR)/include -DBL_SPACEDIM=3 -DAMREX_SPACEDIM=3 -O2
LDFLAGS= -L. -L$(AMREX_DIR)/lib -lamrex -lmpichf90 -lgfortran -lstdc++ -std=c++11 -lm -Wl,-flat_namespace -Wl,-commons,use_dylibs -lmpifort -lmpi -lpmpi  -lgfortran -lquadmath
