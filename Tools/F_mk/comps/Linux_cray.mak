    CC  := cc -target=linux
    FC  := ftn -target=linux
    F90 := ftn -target=linux

    FFLAGS   += -J $(mdir) -I $(mdir) -em
    F90FLAGS += -J $(mdir) -I $(mdir) -em

    ifdef NDEBUG
      FFLAGS   += -O2
      F90FLAGS += -O2
    else
      FFLAGS   += -g -O0
      F90FLAGS += -g -O0
    endif

    ifndef OMP
      FFLAGS   += -h noomp
      F90FLAGS += -h noomp
    endif
