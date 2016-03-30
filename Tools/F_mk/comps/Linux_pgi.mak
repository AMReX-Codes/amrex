
    CC  := pgcc
    CXX := pgCC
    FC  := pgf95
    F90 := pgf95
    FFLAGS   += -module $(mdir) -I$(mdir) 
    F90FLAGS += -module $(mdir) -I$(mdir)

    FCOMP_VERSION := $(shell $(FC) -V 2>&1 | grep 'target')

    ifdef OMP
      F90FLAGS += -mp=nonuma -Minfo=mp
      FFLAGS += -mp=nonuma -Minfo=mp
      CFLAGS += -mp=nonuma -Minfo=mp
      CXXFLAGS += -mp=nonuma -Minfo=mp
    endif
    
    ifdef ACC
      F90FLAGS += -acc -Minfo=acc
      FFLAGS += -acc -Minfo=acc
      CFLAGS += -acc -Minfo=acc
      CXXFLAGS += -acc -Minfo=acc
    else
      F90FLAGS += -noacc
      FFLAGS += -noacc
      CFLAGS += -noacc
      CXXFLAGS += -noacc
    endif

    ifdef NDEBUG
      FFLAGS   += -O2
      F90FLAGS += -O2
      CFLAGS   += -O2
      CXXFLAGS += -O2
    else
      FFLAGS   += -g
      F90FLAGS += -g
      CFLAGS   += -g
      CXXFLAGS += -g
    endif

    xtr_libraries += -pgcpplibs
