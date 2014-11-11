#
# Makefile for DIM dimensional BoxLib wrapper.
#

#
# Configuration
#

TOP            = ../../..
BOXLIB_HOME    = ../..

PRECISION      = DOUBLE
DEBUG	       = FALSE
PROFILE        = FALSE
USE_MPI        = TRUE
USE_OMP        = FALSE

NEEDSCHEM      = FALSE
NEEDSGEOM      = FALSE

ifeq ($(NEEDSCHEM), TRUE)
  COMBUSTION_HOME = $(TOP)/Combustion
  CHEMISTRY_DIR   = $(COMBUSTION_HOME)/Chemistry

  CHEMISTRY_MODEL=GRI12
  CHEMISTRY_MODEL=INERT30
  CHEMISTRY_MODEL=CH4-2STEP
  CHEMISTRY_MODEL=PROPANE
  CHEMISTRY_MODEL=CHEMH
  CHEMISTRY_MODEL=LIDRYMOD
  CHEMISTRY_MODEL=DRM19
  CHEMISTRY_MODEL=GRI30_noArN
  CHEMISTRY_MODEL=LUDME
endif

#
# Definitions
#

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

FFLAGS   += -fPIC
fFLAGS   += -fPIC
CFLAGS   += -fPIC
CXXFLAGS += -fPIC

#
# Paths
#

Bdirs   += $(BOXLIB_HOME)/Src/C_BaseLib
ifeq ($(NEEDSCHEM), TRUE)
  Bdirs += $(COMBUSTION_HOME)/Chemistry/src
endif
ifeq (${NEEDSGEOM}, TRUE)
  Bdirs   += ${BOXLIB_HOME}/Src/C_BoundaryLib
endif

Bpack	+= $(foreach dir, $(Bdirs), $(dir)/Make.package)
Blocs	+= $(foreach dir, $(Bdirs), $(dir))

include $(Bpack)

INCLUDE_LOCATIONS += src $(Blocs)
VPATH_LOCATIONS   += src $(Blocs)

ifeq (${NEEDSCHEM}, TRUE)
  ifeq (${CHEMISTRY_MODEL}, DRM19)
    CHEM_MECHFILE = drm19.c
    VPATH_LOCATIONS += ${CHEMISTRY_DIR}/data/gri/PMFs
    VPATH_LOCATIONS += ${CHEMISTRY_DIR}/data/gri
  endif
  ifeq (${CHEMISTRY_MODEL}, CHEMH)
    CHEM_MECHFILE = chem-H.c
    VPATH_LOCATIONS += ${CHEMISTRY_DIR}/data/chem-H/PMFs
    VPATH_LOCATIONS += ${CHEMISTRY_DIR}/data/chem-H
  endif
  ifeq (${CHEMISTRY_MODEL}, LIDRY)
    CHEM_MECHFILE = LiDryer.c
    VPATH_LOCATIONS += ${CHEMISTRY_DIR}/data/LiDryer/PMFs
    VPATH_LOCATIONS += ${CHEMISTRY_DIR}/data/LiDryer
  endif
  ifeq (${CHEMISTRY_MODEL}, LIDRYMOD)
    CHEM_MECHFILE = LiDryerMOD.c
    VPATH_LOCATIONS += ${CHEMISTRY_DIR}/data/LiDryer/PMFs
    VPATH_LOCATIONS += ${CHEMISTRY_DIR}/data/LiDryer
  endif
  ifeq (${CHEMISTRY_MODEL}, GRI30_noArN)
    CHEM_MECHFILE = grimech30-noArN.c
    VPATH_LOCATIONS += ${CHEMISTRY_DIR}/data/gri/PMFs
    VPATH_LOCATIONS += ${CHEMISTRY_DIR}/data/gri
  endif
  ifeq (${CHEMISTRY_MODEL}, GLARSKEL)
    CHEM_MECHFILE = glarSkel.c
    VPATH_LOCATIONS += ${CHEMISTRY_DIR}/data/glar/PMFs
    VPATH_LOCATIONS += ${CHEMISTRY_DIR}/data/glar
  endif
  ifeq (${CHEMISTRY_MODEL}, LUDME)
    CHEM_MECHFILE = LuDME.c
    VPATH_LOCATIONS += ${CHEMISTRY_DIR}/data/Lu/PMFs
    VPATH_LOCATIONS += ${CHEMISTRY_DIR}/data/Lu
  endif
  cEXE_sources += ${CHEM_MECHFILE}
endif

vpath %.c   $(VPATH_LOCATIONS)
vpath %.cpp $(VPATH_LOCATIONS)
vpath %.h   $(VPATH_LOCATIONS)
vpath %.H   $(VPATH_LOCATIONS)
vpath %.F   $(VPATH_LOCATIONS)
vpath %.f   $(VPATH_LOCATIONS)
vpath %.f90 $(VPATH_LOCATIONS)

#
# Python and NumPy
#

PYINCLUDE := $(shell python -c 'import distutils.sysconfig; print distutils.sysconfig.get_python_inc()')
NPINCLUDE := $(shell python -c 'import numpy; print numpy.get_include()')
PYLIBS    := $(shell python-config --ldflags)

INCLUDE_LOCATIONS += $(PYINCLUDE) $(NPINCLUDE)

FORTLIBS =
ifeq ($(FCOMP), gfortran)
  __gcc_lib_dir := $(dir $(shell gfortran -print-file-name=libgfortran.a))
 FORTLIBS += -L$(__gcc_lib_dir) -lgfortran
endif

#
# Rules
#

WRAPPER = src/boxlib_wrap_$(DIM).cpp
PYSO    = $(OUT)/_bl$(DIM).so

all: $(PYSO)

$(PYSO): $(objForExecs) $(objEXETempDir)/boxlib_wrap_$(DIM).o
	$(CXX) -shared -o $@ $^ ${SHARED_LIBRARIES} ${PYLIBS} ${FORTLIBS}

wrapper: $(WRAPPER)
$(WRAPPER): swig/boxlib.i
	swig -DDIM$(DIM) -python -c++ -Iswig -Icsrc $(includes) -o $@ -outdir boxlib $<

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
