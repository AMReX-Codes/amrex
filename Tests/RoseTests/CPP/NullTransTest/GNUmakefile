
BOXLIB_HOME ?= $(HOME)/Development/BoxLib

BL_NOFAST=TRUE

EBASE = RoseNullTest


DEBUG	  = FALSE
DEBUG	  = TRUE
PROFILE   = TRUE
DIM       = 3
USE_MPI   = FALSE
USE_OMP = FALSE

USE_PARTICLES = TRUE

#COMP      = g++
#FCOMP     = gfortran

COMP      = rose
FCOMP     = rose
ROSECOMP  = identityTranslator
ROSEFCOMP = identityTranslator

HERE = .

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package
include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
include $(BOXLIB_HOME)/Src/C_AMRLib/Make.package


fFLAGS += -ffixed-line-length-none
FFLAGS += -ffixed-line-length-none
ROSECXXFLAGS +=

DEFINES += -DBL_PARALLEL_IO

CEXE_sources += RoseNullTest.cpp
FEXE_sources += PROB_$(DIM)D.F

#INCLUDE_LOCATIONS += $(HERE)
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Extern/amrdata
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BoundaryLib

PATHDIRS  = $(HERE)
PATHDIRS += $(BOXLIB_HOME)/Src/C_BaseLib
PATHDIRS += $(BOXLIB_HOME)/Src/C_AMRLib
PATHDIRS += $(BOXLIB_HOME)/Src/Extern/amrdata
PATHDIRS += $(BOXLIB_HOME)/Src/C_BoundaryLib

vpath %.c   $(PATHDIRS)
vpath %.h   $(PATHDIRS)
vpath %.f   $(PATHDIRS)
vpath %.F   $(PATHDIRS)
vpath %.h   $(PATHDIRS)
vpath %.H   $(PATHDIRS)
vpath %.cpp $(PATHDIRS)
vpath %.f90 $(PATHDIRS)

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules

