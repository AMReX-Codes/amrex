

    ifeq ($(findstring titan, $(HOST)), titan)
        #On Crays like Titan, you need Cray wrappers even for non-Cray compiler
        CXX := CC
        CC  := cc
        FC  := ftn
        F90 := ftn
    else
        CC  := pgcc
        CXX := pgc++
        FC  := pgf95
        F90 := pgf95
    endif
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
      # Disable debug symbols on PGI for now
      ifndef ACC
        FFLAGS   += -gopt -O2
        F90FLAGS += -gopt -O2
        CFLAGS   += -gopt -O2
        CXXFLAGS += -gopt -O2
      endif
    else
      FFLAGS   += -g
      F90FLAGS += -g
      CFLAGS   += -g
      CXXFLAGS += -g
    endif

    CXXFLAGS += --c++11

    ifneq ($(findstring titan, $(HOST)), titan)
        #The wrappers should pick this up on Titan, so don't add it in that case.
        LDFLAGS += -pgc++libs
    endif
