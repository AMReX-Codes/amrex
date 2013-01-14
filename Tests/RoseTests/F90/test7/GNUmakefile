# To compile without ROSE, type "make ROSE= "

F90 = gfortran
F90FLAGS = 

ROSECOMP = identityTranslator
ROSE_F90FLAGS = -rose:verbose 1

ROSE := t

ifndef ROSE
   COMPILE.f90 = $(F90) $(F90FLAGS) -c
else
   COMPILE.f90 = $(ROSECOMP) $(F90FLAGS) $(ROSE_F90FLAGS) -c
endif

LINK.f90 = $(F90) $(F90FLAGS)

default: test7.exe
	@echo SUCCESS

test7.exe: main.o probin.o runtime.o
	$(LINK.f90) -o $@ $^

main.o: main.f90 probin.o
	$(COMPILE.f90) -o $@ $<

probin.o: probin.f90
	$(COMPILE.f90) -o $@ $<

runtime.o: runtime.f90 probin.o
	$(COMPILE.f90) -o $@ $<

clean::
	$(RM) *.o *.mod *.rmod rose_*.f90

realclean:: clean
	$(RM) test7.exe
