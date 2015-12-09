# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME ?= ../..

NDEBUG    := t
MPI       := t
OMP       :=
PROF      :=
COMP      := gfortran
MKVERBOSE := t

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

include ./GPackage.mak
include ../../Source/GPackage.mak

VPATH_LOCATIONS += .	../../Source

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

main.$(suf).exe: $(objects) 
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)
	@echo SUCCESS

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak
