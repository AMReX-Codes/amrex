# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME ?= $(HOME)/Development/BoxLib/

NDEBUG    := t
MPI       := t
OMP       := t
PROF      :=
COMP      := gfortran
MKVERBOSE := t

#USE_IPM   := t

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

include ./GPackage.mak
VPATH_LOCATIONS += .

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

#libraries += $(IPM_GNU)

main.$(suf).exe: $(objects) 
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak
