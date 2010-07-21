  F_C_LINK := LOWERCASE
  COMP = xlf
  ifdef OMP
    rsuf := _r
  endif
  F90 = xlf95$(rsuf)
  FC  = xlf95$(rsuf)
  F90FLAGS += -qnosave -qmoddir=$(mdir) -I$(mdir) -qsuffix=f=f90
  FFLAGS   += -qnosave -qmoddir=$(mdir) -I$(mdir) -qsuffix=f=f -qfixed=72
  ifdef NDEBUG
    ifdef GPROF
      FFLAGS   += -O2
      F90FLAGS += -O2
    else
      FFLAGS   += -O3 -qstrict -qtune=auto -qarch=auto -qcache=auto -NS5000
      F90FLAGS += -O3 -qstrict -qtune=auto -qarch=auto -qcache=auto -NS5000
    endif
  else
    FFLAGS += -g
    FFLAGS += -C
    FFLAGS += -qinitauto=FF 
    FFLAGS += -qlanglvl=95std
    F90FLAGS += -g
    F90FLAGS += -C
    F90FLAGS += -qinitauto=FF
    F90FLAGS += -qlanglvl=95std
  endif
  #
  # You might need the following on seaborg:
  #
  # LDFLAGS += -bmaxdata:0x80000000

  ifdef OMP
    FFLAGS += -qsmp=omp
    F90FLAGS += -qsmp=omp
  endif
  ifdef GPROF
    FFLAGS += -pg
    F90FLAGS += -pg
  endif
