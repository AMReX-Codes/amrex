  COMP = cftn
  F90 := ftn
  FC  := ftn
  CC   := cc
  FFLAGS   =
  F90FLAGS =
  FFLAGS   += -p $(mdir)  -J $(mdir) -e m
  F90FLAGS += -p $(mdir)  -J $(mdir) -e m
  ifndef NDEBUG
    FFLAGS   += -g
    F90FLAGS += -g
  endif
  f2kcli_suf := _crayx1
