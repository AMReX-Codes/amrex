

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
      F90FLAGS += -acc -Minfo=acc -ta=tesla:cuda$(CUDA_VERSION)
      FFLAGS += -acc -Minfo=acc -ta=tesla:cuda$(CUDA_VERSION)
      CFLAGS += -acc -Minfo=acc -ta=tesla:cuda$(CUDA_VERSION)
      CXXFLAGS += -acc -Minfo=acc -ta=tesla:cuda$(CUDA_VERSION)
    else
      F90FLAGS += -noacc
      FFLAGS += -noacc
      CFLAGS += -noacc
      CXXFLAGS += -noacc
    endif

    ifdef CUDA
      F90FLAGS += -Mcuda=cuda$(CUDA_VERSION)
      FFLAGS += -Mcuda=cuda$(CUDA_VERSION)
      CFLAGS += -Mcuda=cuda$(CUDA_VERSION)
      CXXFLAGS += -Mcuda=cuda$(CUDA_VERSION)
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

    ifeq ($(findstring summitdev, $(HOST)), summitdev)
       libraries += -L /sw/summitdev/gcc/5.4.0new/lib64/ -latomic
    endif
