    CC  := cc -target=linux
    FC  := ftn -target=linux
    F90 := ftn -target=linux

    FFLAGS   += -J $(mdir) -I $(mdir) -em
    F90FLAGS += -J $(mdir) -I $(mdir) -em

    ifdef NDEBUG
      FFLAGS   += -O 2,ipa4
      F90FLAGS += -O 2,ipa4
    else
      FFLAGS   += -g
      F90FLAGS += -g
    endif

    ifndef OMP
      FFLAGS   += -h noomp
      F90FLAGS += -h noomp
    endif
