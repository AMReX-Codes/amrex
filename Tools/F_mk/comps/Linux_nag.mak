    FC  = nf95
    F90 = nf95
    CC  = gcc
    F90FLAGS += -mdir $(mdir) -I $(mdir)
    FFLAGS   += -mdir $(mdir) -I $(mdir)
    FFLAGS   += -w=x77 -fixed
    CFLAGS += -Wall
#   F90FLAGS += -Oassumed=always_contig
    f2kcli_suf := _nag
    ifdef NDEBUG
      FFLAGS += -O4
      F90FLAGS += -O4
    else
      FFLAGS += -C=all
      F90FLAGS += -C=all
      FFLAGS += -g
      F90FLAGS += -g
      FFLAGS += -nan
      F90FLAGS += -nan
      FFLAGS += -gline
      F90FLAGS += -gline
    endif
    ifdef GPROF
      FFLAGS += -pg
      F90FLAGS += -pg
      CFLAGS += -pg
    endif
