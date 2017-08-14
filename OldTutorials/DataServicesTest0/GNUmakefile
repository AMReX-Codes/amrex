AMREX_HOME = ../../../amrex

EBASE = DataServicesTest0

DEBUG	  = FALSE
DEBUG	  = TRUE

PROFILE = FALSE
#PROFILE = TRUE
#TRACE_PROFILE = TRUE
#COMM_PROFILE = TRUE 
#MEM_PROFILE = TRUE 

DIM       = 3
USE_MPI   = FALSE
USE_MPI   = TRUE
COMP      = gnu

HERE = .

#DEFINES += -DBL_NOLINEVALUES

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += $(EBASE).cpp
include $(AMREX_HOME)/Src/Extern/ProfParser/Make.package
include $(AMREX_HOME)/Src/Base/Make.package

INCLUDE_LOCATIONS += $(HERE)
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/ProfParser
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata

include $(AMREX_HOME)/Src/Extern/amrdata/Make.package

PATHDIRS  = $(HERE)
PATHDIRS += $(AMREX_HOME)/Src/Base
PATHDIRS += $(AMREX_HOME)/Src/Extern/amrdata
PATHDIRS += $(AMREX_HOME)/Src/Extern/ProfParser


vpath %.f   $(PATHDIRS)
vpath %.F   $(PATHDIRS)
vpath %.f90 $(PATHDIRS)
vpath %.h   $(PATHDIRS)
vpath %.H   $(PATHDIRS)
vpath %.cpp $(PATHDIRS)

all: $(executable)

$(executable):

BLProfParser.tab.H BLProfParser.tab.cpp: $(AMREX_HOME)/Src/Extern/ProfParser/BLProfParser.y
	cat $(AMREX_HOME)/Src/Extern/ProfParser/BLProfParser.y $(SED0) $(SED1) > BLProfParserNC.y
	bison --defines=BLProfParser.tab.H --output=BLProfParser.tab.cpp \
		BLProfParserNC.y
	rm BLProfParserNC.y


BLProfParser.lex.yy.cpp: BLProfParser.tab.H $(AMREX_HOME)/Src/Extern/ProfParser/BLProfParser.l
	flex --outfile=BLProfParser.lex.yy.cpp $(AMREX_HOME)/Src/Extern/ProfParser/BLProfParser.l


cclean:
	rm -f amrprofparser BLProfParser.lex.yy.cpp BLProfParser.tab.cpp \
              BLProfParser.tab.H core

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

