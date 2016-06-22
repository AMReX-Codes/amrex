NDEBUG := t
MPI    :=
OMP    :=
MKVERBOSE :=t 

COMP := gfortran


# include the main Makefile stuff
include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

# to make python libraries, we need to do -fPIC (in gfortran)
F90FLAGS += -fPIC
FFLAGS += -fPIC
CFLAGS += -fPIC
CXXFLAGS += -fPIC

# core BoxLib directories
BOXLIB_CORE := Src/F_BaseLib

Fmpack := $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir))


# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files
VPATH_LOCATIONS += $(Fmlocs)

all: python_module

python_module: $(objects)
	f2py --fcompiler=$(F90) --f90flags="-J $(mdir)" -c fsnapshot.f90 -m fsnapshot $(objects) -lstdc++

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak


clean::
	$(RM) fsnapshot.so


# for debugging.  To see the value of a Makefile variable,                      
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will                        
# print out the value.                                                          
print-%: ; @echo $* is $($*) 
