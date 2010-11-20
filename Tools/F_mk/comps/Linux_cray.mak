    CC  := cc -target=linux
    FC  := ftn -target=linux
    F90 := ftn -target=linux

    FFLAGS   += -J $(mdir) -I $(mdir) -em
    F90FLAGS += -J $(mdir) -I $(mdir) -em

    ifdef NDEBUG
      FFLAGS   += -O 1
      F90FLAGS += -O 1
      CFLAGS   += -O 1
    else
      FFLAGS   += -g -O0
      F90FLAGS += -g -O0
      CFLAGS   += -g -O0
    endif

    ifndef OMP
      FFLAGS   += -h noomp
      F90FLAGS += -h noomp
    endif
