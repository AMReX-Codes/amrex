# To compile without ROSE, type "make ROSE= "

F90 = gfortran
F90FLAGS = -fopenmp

ROSECOMP = identityTranslator
ROSE_F90FLAGS = -rose:verbose 1 -rose:openmp

ROSE := t

ifndef ROSE
   COMPILE.f90 = $(F90) $(F90FLAGS) -c
else
   COMPILE.f90 = $(ROSECOMP) $(F90FLAGS) $(ROSE_F90FLAGS) -c
endif

LINK.f90 = $(F90) $(F90FLAGS)

default: test6.exe
	@echo SUCCESS

test6.exe: main.o omp.o
	$(LINK.f90) -o $@ $^

main.o: main.f90 omp.o
	$(COMPILE.f90) -o $@ $<

omp.o: omp.f90
	$(COMPILE.f90) -o $@ $<

clean::
	$(RM) *.o *.mod *.rmod rose_*.f90

realclean:: clean
	$(RM) test6.exe
