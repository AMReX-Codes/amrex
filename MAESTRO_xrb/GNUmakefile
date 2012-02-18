NDEBUG := t
MPI    :=
OMP    :=
MKVERBOSE :=t 
COMP := gfortran

EOS := extern/EOS/helmeos
NETWORK := extern/networks/null

# some routines need an eos/network (i.e. to compute thermodynamic 
# quantities.  If that is the case, set NEED_EOS_NETWORK := t
NEED_EOS_NETWORK := t

# define the location of the fParallel root directory
FPARALLEL := ../..


# include the main Makefile stuff
include $(FPARALLEL)/mk/GMakedefs.mak

# define the packages to build these tools
Fmdirs := boxlib \
          extern/constants \
	  extern/VODE \
	  extern/BLAS \
	  $(NETWORK) \
	  $(EOS) 

# directories containing files that are 'include'-d via Fortran
Fmincludes := $(EOS)

# ifdef NEED_EOS_NETWORK
#   Fmdirs += extern/EOS/helmeos \
#             extern/networks/triple_alpha \
#             extern/VODE \
# 	    extern/BLAS

#   Fmincludes += extern/EOS/helmeos
# endif

Fmpack := $(foreach dir, $(Fmdirs), $(FPARALLEL)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(Fmdirs), $(FPARALLEL)/$(dir))
Fmincs := $(foreach dir, $(Fmincludes), $(FPARALLEL)/$(dir))

# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files
VPATH_LOCATIONS += $(Fmlocs)

# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)

#programs += fbuoyancy
programs += fad_excess
#programs += fspeciesmass

all: $(pnames)

include $(FPARALLEL)/mk/GMakerules.mak

%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else	
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif



