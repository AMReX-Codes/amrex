# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME ?= $(HOME)/Development/BoxLib/

NDEBUG    := t
MPI       := 
OMP       := 
PROF      := t
COMP      := gfortran
MKVERBOSE := t

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

f90sources += advance.f90
f90sources += main.f90


ifdef OMP
  override CXXOPTF  = -O3 -msse3 -fopenmp
  override COPTF    = -O3 -msse3 -fopenmp
  override FOPTF    = -J$(mdir) -I $(mdir) -O3 -msse3 -fno-range-check -fopenmp
  override fOPTF    = -J$(mdir) -I $(mdir) -O3 -msse3 -fno-range-check -fopenmp
  override F90FLAGS = -J$(mdir) -I $(mdir) -O3 -msse3 -fno-range-check -fopenmp
else
  override CXXOPTF  = -O3 -msse3
  override COPTF    = -O3 -msse3
  override FOPTF    = -J$(mdir) -I $(mdir) -O3 -msse3 -fno-range-check
  override fOPTF    = -J$(mdir) -I $(mdir) -O3 -msse3 -fno-range-check
  override F90FLAGS = -J$(mdir) -I $(mdir) -O3 -msse3 -fno-range-check
endif


VPATH_LOCATIONS += .

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

main.$(suf).exe: $(objects) 
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak
