NDEBUG :=
MPI    :=
OMP    :=

COMP := g95
COMP := Lahey
COMP := Intel

FPARALLEL := ../..

include $(FPARALLEL)/mk/GMakedefs.mak

f90sources += main.f90
f90sources += t_main.f90

include $(FPARALLEL)/boxlib/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/boxlib

main.$(suf).exe: $(objects)
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(FPARALLEL)/mk/GMakerules.mak

