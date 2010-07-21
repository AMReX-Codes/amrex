  ifdef OMP
    FC  := xlf95_r
    F90 := xlf95_r
    CC  := xlc_r

    F90FLAGS += -qsmp=noauto:omp
    FFLAGS   += -qsmp=noauto:omp
    CFLAGS   += -qsmp=noauto:omp
  else
    FC  := xlf95
    F90 := xlf95
    CC  := xlc
  endif 

  F90FLAGS += -I $(mdir) -qmoddir=$(mdir)
  FFLAGS   += -I $(mdir) -qmoddir=$(mdir)
  CFLAGS += 

  F_C_LINK := LOWERCASE
  ifdef NDEBUG
    F90FLAGS += -O 
    FFLAGS += -O 
    CFLAGS += -O
  else
    F90FLAGS += -g 
    FFLAGS += -g 
    CFLAGS += -g
  endif
