AMREX_HOME ?= ../../../amrex

EBASE = DataServicesTest0

DEBUG	  = FALSE
DEBUG	  = TRUE
PROFILE   = FALSE
DIM       = 3
USE_MPI     =FALSE
USE_MPI     =TRUE
COMP      = gnu

HERE = .

DEFINES += -DBL_NOLINEVALUES

include $(AMREX_HOME)/Tools/GNUMake/Make.defs


CEXE_sources += $(EBASE).cpp

include $(AMREX_HOME)/Src/Base/Make.package

INCLUDE_LOCATIONS += $(HERE)
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata

include $(AMREX_HOME)/Src/Extern/amrdata/Make.package

PATHDIRS  = $(HERE)
PATHDIRS += $(AMREX_HOME)/Src/Base
PATHDIRS += $(AMREX_HOME)/Src/Extern/amrdata


vpath %.f   $(PATHDIRS)
vpath %.F   $(PATHDIRS)
vpath %.f90 $(PATHDIRS)
vpath %.h   $(PATHDIRS)
vpath %.H   $(PATHDIRS)
vpath %.cpp $(PATHDIRS)

all: $(executable)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

