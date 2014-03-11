PYBOXLIB ?= $(BOXLIB_HOME)/Src/Python

PYINCLUDE := -I$(shell python -c 'import distutils.sysconfig; print distutils.sysconfig.get_python_inc()')
NPINCLUDE := -I$(shell python -c 'import numpy; print numpy.get_include()')/numpy

CFLAGS += $(PYINCLUDE) $(NPINCLUDE)

csources   += boxlib_numpy_c.c
f90sources += blobjects.f90
f90sources += fboxlib.f90
f90sources += boxlib_numpy_f.f90

PYFBOXLIB = libpyfboxlib.so

VPATH_LOCATIONS += $(PYBOXLIB)/fsrc
