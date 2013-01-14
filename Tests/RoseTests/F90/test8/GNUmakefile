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

default: test8.exe
	@echo SUCCESS

test8.exe: main.o init_data.o
	$(LINK.f90) -o $@ $^

main.o: main.f90 init_data.o
	$(COMPILE.f90) -o $@ $<

init_data.o: init_data.f90
	$(COMPILE.f90) -o $@ $<

clean::
	$(RM) *.o *.mod *.rmod rose_*.f90

realclean:: clean
	$(RM) test8.exe
