# AMREX_HOME defines the directory in which we will find all the BoxLib code
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= $(HOME)/Development/amrex/

NDEBUG    := t
MPI       := t
OMP       := t
PROF      :=
COMP      := gfortran
MKVERBOSE := t

#USE_IPM   := t

include $(AMREX_HOME)/Tools/F_mk/GMakedefs.mak

include ./GPackage.mak
VPATH_LOCATIONS += .

include $(AMREX_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(AMREX_HOME)/Src/F_BaseLib

#libraries += $(IPM_GNU)

main.$(suf).exe: $(objects) 
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(AMREX_HOME)/Tools/F_mk/GMakerules.mak
