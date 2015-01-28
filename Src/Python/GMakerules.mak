F90FLAGS += -fPIC
CFLAGS   += -fPIC

NPINCLUDE:=-I$(shell python -c 'import distutils.sysconfig; print distutils.sysconfig.get_python_inc()')
NPINCLUDE+=-I$(shell python -c 'import numpy; print numpy.get_include()')/numpy

PYLDFLAGS = $(shell python-config --ldflags)
$(PYFBOXLIB): $(objects)
	$(F90) $(F90FLAGS) -shared -o $@ $^ $(PYLDFLAGS)

$(PYCBOXLIB): $(PYBOXLIB)/fsrc/boxlib_numpy_c.c
	$(CC) $(CFLAGS) $(NPINCLUDE) -shared -o $@ $^ -L. -lpyfboxlib $(PYLDFLAGS)
