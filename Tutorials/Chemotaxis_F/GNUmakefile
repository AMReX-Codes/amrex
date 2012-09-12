# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME ?= $(HOME)/Development/BoxLib/

NDEBUG    :=
MPI       :=
OMP       :=
PROF      :=
COMP      := gfortran
MKVERBOSE := 

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

include ./GPackage.mak
VPATH_LOCATIONS += .

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

#include $(BOXLIB_HOME)/Src/Python/GPackage.mak
#all: $(PYBOXLIBS)

main.$(suf).exe: $(objects) 
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak
#include $(BOXLIB_HOME)/Src/Python/GMakerules.mak
