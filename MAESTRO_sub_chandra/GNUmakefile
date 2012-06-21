NDEBUG := t
MPI    :=
OMP    :=

MKVERBOSE :=t 

COMP := gfortran

# some routines need an eos/network (i.e. to compute thermodynamic 
# quantities.  If that is the case, set NEED_EOS_NETWORK := t
NEED_EOS_NETWORK := 

# define the location of the MAESTRO root directory
MAESTRO_TOP_DIR := $(MAESTRO_HOME)


# include the main Makefile stuff
include  $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

# core BoxLib directories
BOXLIB_CORE := Src/F_BaseLib

# other packages needed for data_processing
Fmdirs :=


# directories containing files that are 'include'-d via Fortran
Fmincludes :=

ifdef NEED_EOS_NETWORK
  Fmdirs += Microphysics/EOS/helmeos \
            Microphysics/networks/ignition_simple \
            Util/VODE

  Fmincludes += Microphysics/helmeos
endif

Fmpack := $(foreach dir, $(Fmdirs), $(MAESTRO_TOP_DIR)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(Fmdirs), $(MAESTRO_TOP_DIR)/$(dir))
Fmincs := $(foreach dir, $(Fmincludes), $(MAESTRO_TOP_DIR)/$(dir))

Fmpack += $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir))


# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files
VPATH_LOCATIONS += $(Fmlocs)

# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)


programs += fsubchandra


all: $(pnames)

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak


%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else	
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif



