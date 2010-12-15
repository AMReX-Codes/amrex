  FC  := gfortran
  F90 := gfortran
  CC  := gcc

  F90FLAGS += -J$(mdir) -I $(mdir)
  FFLAGS   += -J$(mdir) -I $(mdir)
  CFLAGS   += -Wall

  ifdef NDEBUG
    F90FLAGS += -O -fno-range-check 
    FFLAGS   += -O -fno-range-check 
    CFLAGS   += -O 
  else
    F90FLAGS += -g -fno-range-check -O -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid -finit-real=nan
    FFLAGS   += -g -fno-range-check -O -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid -finit-real=nan
    CFLAGS   += -g -O
  endif

  ifdef OMP
    F90FLAGS += -fopenmp
    FFLAGS   += -fopenmp
    CFLAGS   += -fopenmp
  endif
