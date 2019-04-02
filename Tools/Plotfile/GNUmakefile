AMREX_HOME ?= ../..

DEBUG = FALSE
DIM = 3
COMP = gcc
PRECISION = DOUBLE
USE_MPI = FALSE
USE_OMP = FALSE
USE_CUDA = FALSE

programs ?=

ifeq ($(strip $(programs)),)
  programs += fboxinfo
  programs += fcompare
  programs += fextract
  programs += fextrema
  programs += fnan
  programs += fsnapshot
  programs += ftime
  programs += fvarnames
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

multiple_executables = $(addsuffix .$(machineSuffix).ex, $(programs))
default: $(multiple_executables)

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

clean::
	$(SILENT) $(RM) $(multiple_executables)
