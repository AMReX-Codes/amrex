    CC  := cc -target=linux
    ifdef MPI
      FC  := ftn -target=linux -module $(mdir) -I$(mdir) 
      F90 := ftn -target=linux -module $(mdir) -I$(mdir) 
    else
      FC  := pgf95 -module $(mdir) -I$(mdir) 
      F90 := pgf95 -module $(mdir) -I$(mdir) 
    endif        
    ifdef NDEBUG
      FFLAGS   += -O
      F90FLAGS += -O
    else
      FFLAGS   += -g -Mbounds
      F90FLAGS += -g -Mbounds
    endif
