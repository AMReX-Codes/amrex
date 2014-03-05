
TOP            = ../CCSE
BOXLIB_HOME    = ${TOP}/BoxLib
COMBUSTION_HOME = ${TOP}/Combustion
CHEMISTRY_DIR = ${COMBUSTION_HOME}/Chemistry

PRECISION      = DOUBLE
DEBUG	       = FALSE
DEBUG	       = TRUE
PROFILE        = FALSE
DIM    	       = 3
DIM    	       = 2
USE_MPI        = FALSE
USE_MPI        = TRUE
USE_OMP        = TRUE
USE_OMP        = FALSE
COMP           = g++
FCOMP          = gfortran
NEEDSFORT      = TRUE
NEEDSCHEM      = TRUE
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

include ${BOXLIB_HOME}/Tools/C_mk/Make.defs

# build to link to shared libs
FFLAGS   += -fPIC
fFLAGS   += -fPIC
CFLAGS   += -fPIC
CXXFLAGS += -fPIC

fincludes=${includes}

# Chemistry
ifeq (${NEEDSCHEM}, TRUE)
  Bdirs += ${COMBUSTION_HOME}/Chemistry/src
endif

# BoxLib
Bdirs   += ${BOXLIB_HOME}/Src/C_BaseLib

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
INCLUDE_LOCATIONS += . $(Blocs)
VPATH_LOCATIONS   += . $(Blocs)

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
endif

VPATH_LOCATIONS += .

cEXE_sources += ${CHEM_MECHFILE}
ifeq (${BL_HAS_FORT}, TRUE)
  FEXE_sources += ${EBASE}_F.F
endif

CEXE_sources += support.cpp chemSupport.cpp
CEXE_headers += support.H chemSupport.H


vpath %.c   $(VPATH_LOCATIONS)
vpath %.cpp $(VPATH_LOCATIONS)
vpath %.h   $(VPATH_LOCATIONS)
vpath %.H   $(VPATH_LOCATIONS)
vpath %.F   $(VPATH_LOCATIONS)
vpath %.f   $(VPATH_LOCATIONS)
vpath %.f90 $(VPATH_LOCATIONS)

ifeq ($(MACHINE), Linux)
  ifeq (${WHICHLINUX}, HEDORAH)
    PYTHON_DIR := /usr
    PYTHON_INCLUDES := ${PYTHON_DIR}/include/python2.6
    ifeq ($(USE_MPI), TRUE)
      SHARED_LIBRARIES += /usr/local/lib/libmpichcxx.so \
                          /usr/local/lib/libmpl.so
    endif
    SHARED_LIBRARIES += /usr/lib/libstdc++.so.6 \
                        /usr/lib/gcc/x86_64-linux-gnu/4.4.3/libgfortran.so \
                        /usr/lib/libm.so
    endif
endif
ifeq (${MACHINE}, Darwin)
  PYTHON_DIR := /usr/local/Cellar/python/2.7.6/Frameworks/Python.framework/Versions/Current
  PYTHON_INCLUDES := ${PYTHON_DIR}/Headers
    ifeq ($(USE_MPI), TRUE)
      SHARED_LIBRARIES += /usr/local/openmpi/lib/libmpi.dylib /usr/local/openmpi/lib/libmpi_cxx.dylib
    endif
  SHARED_LIBRARIES += ${PYTHON_DIR}/lib/libpython2.7.dylib \
                      /usr/local/lib/libgfortran.dylib
  PYTHON_INCLUDES += $(shell python -c 'import numpy; print numpy.get_include()')
endif

INCLUDE_LOCATIONS += ${PYTHON_INCLUDES}

_boxlib.so: all $(objForExecs) $(objEXETempDir)/boxlib_wrap.o
	g++ -shared -o _boxlib.so $(objForExecs) $(objEXETempDir)/boxlib_wrap.o ${SHARED_LIBRARIES}

boxlib_wrap.cpp: boxlib.i
	swig -python -c++ -I. $(includes) -o boxlib_wrap.cpp boxlib.i

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules

realclean::
	rm -rf dat dat1 dat2 dat3 _boxlib.so boxlib_wrap.cpp boxlib.py *.pyc d f o

