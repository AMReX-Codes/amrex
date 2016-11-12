    CC  := cc 
    CXX := CC
    FC  := ftn 
    F90 := ftn 

    FCOMP_VERSION := $(shell ftn -V 2>&1 | grep 'Version')

    FFLAGS   += -J $(mdir) -I $(mdir) -em -hlist=a
    F90FLAGS += -J $(mdir) -I $(mdir) -em -hlist=a

    ifdef NDEBUG
      FFLAGS   += -G2 -O 1
      F90FLAGS += -G2 -O 1
      CFLAGS   += -G2 -O 1
      CXXFLAGS += -G2 -O 1
    else
      FFLAGS   += -g -O0
      F90FLAGS += -g -O0
      CFLAGS   += -g -O0
      CXXFLAGS += -g -O0
    endif

    ifdef ACC
      #These are based on Blue Waters suggestions, might need to edit to be more general
      ifndef NDEBUG
        FFLAGS   += -h msgs
        F90FLAGS += -h msgs
        CFLAGS   += -h pragma=msgs
        CXXFLAGS += -h pragma=msgs
      endif
      FFLAGS   += -h acc -fpic -dynamic -lcudart
      F90FLAGS += -h acc -fpic -dynamic -lcudart
      CFLAGS   += -h pragma=acc
      CXXFLAGS += -h pragma=acc
    else
      FFLAGS   += -h noacc
      F90FLAGS += -h noacc
      CFLAGS   += -h nopragma=acc
      CXXFLAGS += -h nopragma=acc
    endif

    ifndef OMP
      FFLAGS   += -h noomp
      F90FLAGS += -h noomp
      CFLAGS   += -h noomp
      CXXFLAGS += -h noomp
    endif
 
    CXXFLAGS += -hstd=c++11
