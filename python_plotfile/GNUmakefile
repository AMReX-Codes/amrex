NDEBUG := t
MPI    :=
OMP    :=
MKVERBOSE :=t 
COMP := gfortran


# some routines need an eos/network (i.e. to compute thermodynamic 
# quantities.  If that is the case, set NEED_EOS_NETWORK := t
NEED_EOS_NETWORK := 

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
Fmdirs := extern/constants 


# directories containing files that are 'include'-d via Fortran
Fmincludes := 

ifdef NEED_EOS_NETWORK
  Fmdirs += extern/EOS/helmeos \
            extern/networks/ignition \
            extern/VODE

  Fmincludes += extern/helmeos
endif

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


