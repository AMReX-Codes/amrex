PYBOXLIB ?= $(BOXLIB_HOME)/Src/Python

libpyfboxlib.so: $(objects)
	mpif90 -shared -o libpyfboxlib.so $(objects)

NPINCLUDE:=-I$(shell python -c 'import distutils.sysconfig; print distutils.sysconfig.get_python_inc()')
NPINCLUDE+=-I$(shell python -c 'import numpy; print numpy.get_include()')/numpy

libpycboxlib.so: $(PYBOXLIB)/src/boxlib_numpy.c
	gcc -fPIC -shared $(NPINCLUDE) -L. -lpyfboxlib -o libpycboxlib.so $(PYBOXLIB)/src/boxlib_numpy.c
