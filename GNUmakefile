BOXLIB_HOME ?= ../BoxLib
PICSAR_HOME = ../picsar

DEBUG	= FALSE
DEBUG	= TRUE

DIM	= 3

COMP    = gcc
FCOMP   = gfortran

USE_PARTICLES = TRUE

PRECISION = DOUBLE

TINY_PROFILE = TRUE

USE_MPI   = TRUE
USE_OMP   = FALSE

EBASE     = main

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

include ./Make.package

include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package
include $(BOXLIB_HOME)/Src/C_ParticleLib/Make.package

include $(PICSAR_HOME)/src/Make.package

all: $(executable) 
	$(SILENT) $(RM) buildInfo.cpp
	@echo SUCCESS

# job_info support
CEXE_sources += buildInfo.cpp
CEXE_headers += $(BOXLIB_HOME)/Tools/C_scripts/buildInfo.H
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Tools/C_scripts

buildInfo.cpp: 
	$(BOXLIB_HOME)/Tools/C_scripts/makebuildinfo_C.py \
          --boxlib_home "$(BOXLIB_HOME)" \
          --COMP "$(COMP)" --COMP_VERSION "$(COMP_VERSION)" \
          --FCOMP "$(FCOMP)" --FCOMP_VERSION "$(FCOMP_VERSION)" \
          --GIT ". $(BOXLIB_HOME) $(PICSAR_HOME)"

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules

clean::
	$(SILENT) $(RM) buildInfo.cpp
