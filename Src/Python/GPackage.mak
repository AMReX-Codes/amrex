PYBOXLIB ?= $(BOXLIB_HOME)/Src/Python

f90sources += blobjects.f90
f90sources += fboxlib.f90
f90sources += boxlib_numpy.f90

PYBOXLIBS= libpyfboxlib.so libpycboxlib.so
VPATH_LOCATIONS += $(PYBOXLIB)/src
