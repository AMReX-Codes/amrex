  COMP = f90
  F90 = f90
  FC  = f90
  ifdef DEBUG
    FFLAGS += -g   -check bounds
    F90FLAGS += -g  -check bounds
  else
    FFLAGS += -fast -inline all
    F90FLAGS += -fast -inline all
  endif
  ifdef OMP
    FFLAGS += -omp
    F90FLAGS += -omp
    ifdef DEBUG
      FFLAGS += -check omp_bindings
      F90FLAGS += -check omp_bindings
    endif
  endif
