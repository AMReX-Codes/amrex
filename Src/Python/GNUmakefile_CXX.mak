
TOP            = $(HOME)/Development
BOXLIB_HOME    = $(TOP)/BoxLib

PRECISION      = DOUBLE
DEBUG	       = FALSE
PROFILE        = FALSE
USE_MPI        = TRUE
USE_OMP        = FALSE

COMP           = g++
FCOMP          = gfortran

NEEDSCHEM      = FALSE

ifeq ($(NEEDSCHEM), TRUE)
  COMBUSTION_HOME = $(TOP)/Combustion
  CHEMISTRY_DIR   = $(COMBUSTION_HOME)/Chemistry

  #
  # Define reaction mechanism & thermal properties.
  #
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

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

# default to using mpi compilers
CXX := mpic++
CC  := mpicc
FC  := mpif90
fC  := mpif90
F90 := mpif90

FFLAGS   += -fPIC
fFLAGS   += -fPIC
CFLAGS   += -fPIC
CXXFLAGS += -fPIC -D__USE_XOPEN2K8

# Chemistry
ifeq ($(NEEDSCHEM), TRUE)
  Bdirs += $(COMBUSTION_HOME)/Chemistry/src
endif

# BoxLib
Bdirs   += $(BOXLIB_HOME)/Src/C_BaseLib

# Remove for now, add later
#ifeq (${NEEDSGEOM}, TRUE)
#  Bdirs   += ${BOXLIB_HOME}/Src/C_BoundaryLib
#endif
#Bdirs   += ${BOXLIB_HOME}/Src/Extern/amrdata
#DEFINES += -DBL_PARALLEL_IO -DBL_NOLINEVALUES
#FEXE_sources += FILCC_${DIM}D.F

Bpack	+= $(foreach dir, $(Bdirs), $(dir)/Make.package)
Blocs	+= $(foreach dir, $(Bdirs), $(dir))

include $(Bpack)
INCLUDE_LOCATIONS += src $(Blocs)
VPATH_LOCATIONS   += src $(Blocs)

# Remove for now, add later
#INCLUDE_LOCATIONS += ${BOXLIB_HOME}/Src/C_AMRLib

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

ifeq (${BL_HAS_FORT}, TRUE)
  FEXE_sources += ${EBASE}_F.F
endif

# CEXE_sources += support.cpp chemSupport.cpp
# CEXE_headers += support.H chemSupport.H

vpath %.c   $(VPATH_LOCATIONS)
vpath %.cpp $(VPATH_LOCATIONS)
vpath %.h   $(VPATH_LOCATIONS)
vpath %.H   $(VPATH_LOCATIONS)
vpath %.F   $(VPATH_LOCATIONS)
vpath %.f   $(VPATH_LOCATIONS)
vpath %.f90 $(VPATH_LOCATIONS)

PYINCLUDE := $(shell python -c 'import distutils.sysconfig; print distutils.sysconfig.get_python_inc()')
NPINCLUDE := $(shell python -c 'import numpy; print numpy.get_include()')

INCLUDE_LOCATIONS += $(PYINCLUDE) $(NPINCLUDE)

WRAPPER = src/boxlib_wrap_$(DIM).cpp
PYSO    = $(OUT)/_bl$(DIM).so

all: $(PYSO)

#wrapper: $(WRAPPER)

#$(WRAPPER): swig/boxlib.i
#	swig -DDIM$(DIM) -python -c++ -Iswig -Icsrc $(includes) -o $@ -outdir boxlib $<

#mpic++ -print-file-name=libstdc++.a

$(PYSO): $(objForExecs) $(objEXETempDir)/boxlib_wrap_$(DIM).o
	$(F90) -shared -o $@ $^ ${SHARED_LIBRARIES} -lstdc++

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules


