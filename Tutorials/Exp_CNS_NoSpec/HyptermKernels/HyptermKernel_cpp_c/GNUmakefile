BOXLIB_HOME ?= ../../../../../BoxLib

BL_NOFAST=TRUE

#
# Variables for the user to set ...
#
########## define EBASE and BL_HAS_FORT here

EBASE = HyptermOnly

BL_HAS_FORT = TRUE

DEFINES += -DBL_NOBENCHMAIN

DEBUG	  = TRUE
DEBUG	  = FALSE
PROFILE   = FALSE
USE_THREADS =TRUE
USE_THREADS =FALSE
DIM       = 2
DIM       = 3
USE_MPI     =FALSE
USE_MPI     =TRUE
AIXLARGEMEM =TRUE
AIXLARGEMEM =FALSE
COMP      = g++
FCOMP     = gfortran
USE_ARRAYVIEW = FALSE
USE_OMP = TRUE

#DEBUG	  = TRUE
USE_OMP = FALSE
USE_MPI     =FALSE
#FDEBF += -g
#fDEBF += -g
#CDEBF += -g
#CPPDEBF += -g

COPTF += -O3

HERE = .

ifeq ($(EBASE), HyptermOnly)
endif


include $(BOXLIB_HOME)/Tools/C_mk/Make.defs


CEXE_sources += $(EBASE).cpp

include ./Make.package

DEFINES += -DBL_PARALLEL_IO

#FOPTF += -msse3

include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

INCLUDE_LOCATIONS += $(HERE)
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Extern/amrdata
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BoundaryLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Tutorials/Exp_CNS_NoSpec/HyptermKernels/HyptermKernel_cpp_F


include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
INCLUDES += $(BOXLIB_HOME)/Src/C_BoundaryLib

PATHDIRS  = $(HERE)
PATHDIRS += $(BOXLIB_HOME)/Src/C_BaseLib
PATHDIRS += $(BOXLIB_HOME)/Src/C_AMRLib
PATHDIRS += $(BOXLIB_HOME)/Src/Extern/amrdata
PATHDIRS += $(BOXLIB_HOME)/Src/C_BoundaryLib
PATHDIRS += $(BOXLIB_HOME)/Tutorials/Exp_CNS_NoSpec/HyptermKernels/HyptermKernel_c_c
PATHDIRS += $(BOXLIB_HOME)/Tutorials/Exp_CNS_NoSpec/HyptermKernels/HyptermKernel_cpp_F


ifeq (${BL_HAS_FORT}, TRUE)
FEXE_sources += ${EBASE}_F.F
endif

vpath %.f   $(PATHDIRS)
vpath %.F   $(PATHDIRS)
vpath %.h   $(PATHDIRS)
vpath %.H   $(PATHDIRS)
vpath %.cpp $(PATHDIRS)
vpath %.c   $(PATHDIRS)

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules


