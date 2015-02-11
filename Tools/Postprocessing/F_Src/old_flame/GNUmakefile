NDEBUG := t
MPI    :=
OMP    :=
MKVERBOSE :=t 
COMP := Intel

# some routines need an eos/network (i.e. to compute thermodynamic 
# quantities.  If that is the case, set NEED_EOS_NETWORK := t
NEED_EOS_NETWORK := 

# define the location of the fParallel root directory
FPARALLEL := ../..


# include the main Makefile stuff
include $(FPARALLEL)/mk/GMakedefs.mak

# define the packages to build these tools
Fmdirs := boxlib \
          extern/constants

# directories containing files that are 'include'-d via Fortran
Fmincludes := 

ifdef NEED_EOS_NETWORK
  Fmdirs += extern/helmeos \
            extern/networks/ignition \
            extern/VODE

  Fmincludes += extern/helmeos
endif

Fmpack := $(foreach dir, $(Fmdirs), $(FPARALLEL)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(Fmdirs), $(FPARALLEL)/$(dir))
Fmincs := $(foreach dir, $(Fmincludes), $(FPARALLEL)/$(dir))

# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files
VPATH_LOCATIONS += $(Fmlocs)

# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)


#programs += fbubble_position_3d
programs += fbubble_position
#programs += fcusp
#programs += fcylflame
#programs += fflamelength
#programs += finteg
#programs += fturbkin
#programs += fwidth

all: $(pnames)

include $(FPARALLEL)/mk/GMakerules.mak

%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else	
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif



