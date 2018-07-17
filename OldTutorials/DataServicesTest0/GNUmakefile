AMREX_HOME = ../../../amrex

EBASE = DataServicesTest0

DEBUG	  = FALSE
DEBUG	  = TRUE

#TINY_PROFILE=TRUE
#PROFILE = FALSE
#PROFILE = TRUE

#TRACE_PROFILE = TRUE
#COMM_PROFILE = TRUE 
#MEM_PROFILE = TRUE 
USE_PROFPARSER = TRUE

DIM       = 3 
USE_MPI   = FALSE
USE_MPI   = TRUE
COMP      = gnu

HERE = .

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += $(EBASE).cpp

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Extern/amrdata/Make.package

INCLUDE_LOCATIONS += $(HERE)
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata

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

$(executable):

cclean:
	rm -f amrprofparser BLProfParser.lex.yy.cpp BLProfParser.tab.cpp \
              BLProfParser.tab.H core

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

