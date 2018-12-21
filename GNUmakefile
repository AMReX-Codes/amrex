AMREX_INSTALL_DIR = /opt/amrex2d-18.11
DIM = 2
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_ACC = FALSE
COMP = gnu
DEBUG = TRUE
USE_PARTICLES = TRUE
USE_FORTRAN_INTERFACE = TRUE
USE_LINEAR_SOLVERS = TRUE
USE_HYPRE = FALSE
USE_EB = FALSE
AMREX_XSDK = FALSE
ALLOW_DIFFERENT_COMP = FALSE
USE_SENSEI_INSITU = FALSE
TINY_PROFILE = TRUE

AMREX_HOME := $(shell pwd)

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

Pdirs := Base AmrCore Amr Boundary
ifeq ($(USE_PARTICLES),TRUE)
    Pdirs += Particle
endif
ifeq ($(USE_FORTRAN_INTERFACE),TRUE)
    Pdirs += F_Interfaces/Base F_Interfaces/Octree F_Interfaces/AmrCore
endif
ifeq ($(USE_LINEAR_SOLVERS),TRUE)
   Pdirs += LinearSolvers/MLMG LinearSolvers/C_CellMG
   ifeq ($(USE_FORTRAN_INTERFACE),TRUE)
     Pdirs += F_Interfaces/LinearSolvers
   endif
endif
ifeq ($(USE_EB),TRUE)
   Pdirs += EB
endif
ifeq ($(USE_HYPRE),TRUE)
   ifeq ($(USE_LINEAR_SOLVERS),TRUE)
      Pdirs += Extern/HYPRE
   endif
endif
ifeq ($(USE_SENSEI_INSITU),TRUE)
	Pdirs += Extern/SENSEI
endif
Ppack := $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)
include $(Ppack)

all: $(amrexlib)
	@echo SUCCESS

.PHONY: distclean install uninstall

distclean: realclean
	$(SILENT) $(RM) GNUmakefile

install: install_lib install_headers install_fortran_modules

uninstall:
	@echo Uninstalling...
	$(SILENT) $(RM) -r $(AMREX_INSTALL_DIR)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
