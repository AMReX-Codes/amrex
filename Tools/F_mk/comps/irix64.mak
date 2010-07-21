  COMP = f90
  F90  = f90
  FC   = f90
  CC   = cc
  tdir = .
  odir = .
  mdir = .
  ifdef NDEBUG
#   FFLAGS += -O
#   F90FLAGS += -O
      FFLAGS += -Ofast
      F90FLAGS += -Ofast
      LDFLAGS += -Ofast
  else
#   FFLAGS += -C
    FFLAGS += -DEBUG:verbose_runtime=ON
    FFLAGS += -DEBUG:subscript_check=ON
#   FFLAGS += -DEBUG:trap_uninitialized=ON
    FFLAGS += -g
    FFLAGS += -ansi
#   F90FLAGS += -C
    F90FLAGS += -DEBUG:verbose_runtime=ON
    F90FLAGS += -DEBUG:conform_check=ON
    F90FLAGS += -DEBUG:subscript_check=ON
#   F90FLAGS += -DEBUG:trap_uninitialized=ON
    F90FLAGS += -g
    F90FLAGS += -ansi
  endif
  F90FLAGS += -64
  FFLAGS += -64
  CFLAGS += -64
  ifdef OMP
    F90FLAGS += -mp1 
    FFLAGS += -mp1
  endif
