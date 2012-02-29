NDEBUG := t
MPI    :=
OMP    :=
MKVERBOSE :=t 
COMP := gfortran


# define the location of the fParallel root directory
FPARALLEL := ../..


# include the main Makefile stuff
include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

# to make python libraries, we need to do -fPIC (in gfortran)
F90FLAGS += -fPIC
FFLAGS += -fPIC
CFLAGS += -fPIC

# core BoxLib directories
BOXLIB_CORE := Src/F_BaseLib

# other packages needed for data_processing
Fmdirs :=  


# directories containing files that are 'include'-d via Fortran
Fmincludes := 


Fmpack := $(foreach dir, $(Fmdirs), $(FPARALLEL)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(Fmdirs), $(FPARALLEL)/$(dir))
Fmincs := $(foreach dir, $(Fmincludes), $(FPARALLEL)/$(dir))

Fmpack += $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir))



# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files
VPATH_LOCATIONS += $(Fmlocs)

# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)

programs += fsnapshot


all: python_module

python_module: $(objects)
	f2py --fcompiler=gfortran --f90flags="-J t/Linux.gfortran/m/" -c fsnapshot.f90 -m fsnapshot $(objects)

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak


clean::
	$(RM) fsnapshot.so


