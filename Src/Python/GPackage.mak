PYBOXLIB ?= $(BOXLIB_HOME)/Src/Python

PYFBOXLIB = libpyfboxlib.so
PYCBOXLIB = libpycboxlib.so
PYBOXLIBS = $(PYFBOXLIB) $(PYCBOXLIB)

f90sources += blobjects.f90
f90sources += fboxlib.f90
f90sources += boxlib_numpy_f.f90
 
VPATH_LOCATIONS += $(PYBOXLIB)/src
