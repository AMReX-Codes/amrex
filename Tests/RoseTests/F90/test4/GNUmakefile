# To compile without ROSE, type "make ROSE= "

F90 = gfortran
F90FLAGS = 

ROSECOMP = identityTranslator
ROSEFLAGS = -rose:verbose 1

ROSE := t

ifndef ROSE
   COMPILE.f90 = $(F90) $(F90FLAGS) -c
else
   COMPILE.f90 = $(ROSECOMP) $(F90FLAGS) $(ROSEFLAGS) -c
endif

LINK.f90 = $(F90) $(F90FLAGS)

default: main.o
	$(LINK.f90) $^

%.o: %.f90
	$(COMPILE.f90) $^

clean::
	$(RM) *.o *.mod *.rmod rose_*.f90

realclean:: clean
	$(RM) a.out
