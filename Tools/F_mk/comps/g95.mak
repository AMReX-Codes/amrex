  F_C_LINK := DBL_UNDERSCORE
  FC := g95
  F90 := g95
  CC := gcc
# F90FLAGS += -std=f95 -fintrinsic-extensions=iargc,getarg
# FFLAGS   += -std=f95 -fintrinsic-extensions=iargc,getarg  
  F90FLAGS += -fmod=$(mdir) -I $(mdir)
  FFLAGS   += -fmod=$(mdir) -I $(mdir)
# F90FLAGS += -Wall 
# FFLAGS += -Wall 
# CFLAGS += -Wall
# FFLAGS += -ffloat-store
# F90FLAGS += -ffloat-store
  ifdef NDEBUG
    F90FLAGS += -O3 -ffloat-store
    FFLAGS += -O3 -ffloat-store
    CFLAGS += -O3
  else
    F90FLAGS += -g
    F90FLAGS += -fbounds-check
    F90FLAGS += -freal=nan
    FFLAGS += -g
    FFLAGS += -fbounds-check
    FFLAGS += -freal=nan
    CFLAGS += -g
  endif
  ifdef GPROF
    F90FLAGS += -pg
    FFLAGS += -pg
    CFLAGS += -pg
  endif
  ifeq ($(ARCH),Darwin)
    xtr_libraries += -lSystemStubs
  endif
