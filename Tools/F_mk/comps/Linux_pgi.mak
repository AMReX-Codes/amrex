

    ifeq ($(findstring titan, $(HOST)), titan)
        #On Crays like Titan, you need Cray wrappers even for non-Cray compiler
        CXX := CC
        CC  := cc
        FC  := ftn
        F90 := ftn
    else
        CC  := pgcc
        CXX := pgc++
        FC  := pgfortran
        F90 := pgfortran
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
      F90FLAGS += -acc -Minfo=acc -ta=tesla:cuda8.0
      FFLAGS += -acc -Minfo=acc -ta=tesla:cuda8.0
      CFLAGS += -acc -Minfo=acc -ta=tesla:cuda8.0
      CXXFLAGS += -acc -Minfo=acc -ta=tesla:cuda8.0
    else
      F90FLAGS += -noacc
      FFLAGS += -noacc
      CFLAGS += -noacc
      CXXFLAGS += -noacc
    endif

    ifdef CUDA
      F90FLAGS += -Mcuda=cuda8.0
      FFLAGS += -Mcuda=cuda8.0
      CFLAGS += -Mcuda=cuda8.0
      CXXFLAGS += -Mcuda=cuda8.0
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

    ifeq ($(findstring summit, $(HOST)), summit)
       libraries += -L /sw/summitdev/gcc/5.4.0new/lib64/ -latomic
    endif
