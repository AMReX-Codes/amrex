    CC  := cc -target=catamount
    FC  := ftn -target=catamount -module $(mdir) -I$(mdir) 
    F90 := ftn -target=catamount -module $(mdir) -I$(mdir) 
    ifdef NDEBUG
      FFLAGS   += -O
      F90FLAGS += -O
    else
      FFLAGS   += -g
      F90FLAGS += -g
    endif
