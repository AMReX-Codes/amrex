# AMREX_HOME defines the directory in which we will find all the BoxLib code
AMREX_HOME ?= /project/projectdirs/nyx/src/amrex

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ..

EBASE     = AmrDerive

# compilation options
COMP    = gcc

USE_MPI = TRUE
USE_OMP = FALSE

PRECISION = DOUBLE
DEBUG     = TRUE
DEBUG     = FALSE

DIM      = 3

DEFINES += -DBL_NOLINEVALUES

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

Bpack := ./Make.package
Blocs := .

include $(Bpack)
INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

Pdirs   := Base Extern/amrdata
Ppack   += $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)
Plocs   += $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir))

include $(Ppack)
INCLUDE_LOCATIONS += $(Plocs)
VPATH_LOCATIONS   += $(Plocs)

vpath %.c   . $(VPATH_LOCATIONS)
vpath %.cpp . $(VPATH_LOCATIONS)
vpath %.h   . $(VPATH_LOCATIONS)
vpath %.H   . $(VPATH_LOCATIONS)
vpath %.F   . $(VPATH_LOCATIONS)
vpath %.f90 . $(VPATH_LOCATIONS)
vpath %.f   . $(VPATH_LOCATIONS)
vpath %.fi  . $(VPATH_LOCATIONS)

all: $(executable) 
	@echo SUCCESS

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
