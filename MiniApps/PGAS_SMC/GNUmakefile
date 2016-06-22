
BOXLIB_HOME ?= ../..

# compilation options
COMP    = g++
FCOMP   = gfortran

USE_MPI   = TRUE
USE_UPCXX = TRUE

USE_OMP   = FALSE

PROFILE   = FALSE

PRECISION = DOUBLE

DEBUG     = FALSE

DIM      = 3

Bpack := ./Make.package
Blocs := .

EBASE = SMC

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

all: $(executable)
	@echo SUCCESS

include $(Bpack)
INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

Pdirs 	:= C_BaseLib
Ppack	+= $(foreach dir, $(Pdirs), $(BOXLIB_HOME)/Src/$(dir)/Make.package)
Plocs	+= $(foreach dir, $(Pdirs), $(BOXLIB_HOME)/Src/$(dir))

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


include $(BOXLIB_HOME)/Tools/C_mk/Make.rules


#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)

