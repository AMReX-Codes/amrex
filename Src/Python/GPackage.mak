PYBOXLIB ?= $(BOXLIB_HOME)/Src/Python

pybl_sources := $(PYBOXLIB)/src/fboxlib.f90
pybl_pyfs = $(patsubst %.f90,%.pyf,$(addprefix $(tdir)/, $(notdir $(pybl_sources))))

VPATH_LOCATIONS += $(PYBOXLIB)/src
