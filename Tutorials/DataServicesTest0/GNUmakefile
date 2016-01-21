BOXLIB_HOME ?= ../../../BoxLib

EBASE = DataServicesTest0

DEBUG	  = FALSE
DEBUG	  = TRUE
PROFILE   = FALSE
DIM       = 3
USE_MPI     =FALSE
USE_MPI     =TRUE
COMP      = g++
FCOMP     = gfortran

HERE = .

DEFINES += -DBL_NOLINEVALUES

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs


CEXE_sources += $(EBASE).cpp

include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

INCLUDE_LOCATIONS += $(HERE)
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Extern/amrdata

include $(BOXLIB_HOME)/Src/Extern/amrdata/Make.package

PATHDIRS  = $(HERE)
PATHDIRS += $(BOXLIB_HOME)/Src/C_BaseLib
PATHDIRS += $(BOXLIB_HOME)/Src/Extern/amrdata


vpath %.f   $(PATHDIRS)
vpath %.F   $(PATHDIRS)
vpath %.f90 $(PATHDIRS)
vpath %.h   $(PATHDIRS)
vpath %.H   $(PATHDIRS)
vpath %.cpp $(PATHDIRS)

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules

