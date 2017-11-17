    _ifc_version := $(shell ifort -V 2>&1 | grep 'Version')
    ifeq ($(findstring Version 18, $(_ifc_version)), Version 18)
        _comp := Intel18
    else ifeq ($(findstring Version 17, $(_ifc_version)), Version 17)
        _comp := Intel17
    else ifeq ($(findstring Version 16, $(_ifc_version)), Version 16)
        _comp := Intel16
    else ifeq ($(findstring Version 15, $(_ifc_version)), Version 15)
        _comp := Intel15
    else ifeq ($(findstring Version 14, $(_ifc_version)), Version 14)
        _comp := Intel14
    else ifeq ($(findstring Version 13, $(_ifc_version)), Version 13)
        _comp := Intel13
    else
      $(error "$(_ifc_version) of IFC is not supported")
    endif

    F90 := ifort
    FC  := ifort
    CC  := icc
    CXX := icpc

    FCOMP_VERSION := $(shell ifort -V 2>&1 | grep 'Version')

    FFLAGS   += -module $(mdir) -I $(mdir)
    F90FLAGS += -module $(mdir) -I $(mdir) -cxxlib
    CFLAGS   += -std=c99
    CXXFLAGS += -std=c++11

    ifdef MIC
      FFLAGS   += -mmic
      F90FLAGS += -mmic
      CFLAGS   += -mmic
      CXXFLAGS += -mmic
    endif

    ifdef TARGETHOST
      FFLAGS   += -xHost
      F90FLAGS += -xHost
      CFLAGS   += -xHost
      CXXFLAGS += -xHost
    endif

    ifdef OMP
      FFLAGS   += -qopenmp
      F90FLAGS += -qopenmp
      CFLAGS   += -qopenmp
      CXXFLAGS += -qopenmp
    endif

#ifeq ($(GCC_MINOR),$(filter $(GCC_MINOR),4 5))

    ifeq ($(_comp),$(filter $(_comp),Intel17 Intel18))
      ifndef NDEBUG
        F90FLAGS += -g -traceback -O0 -fpe0 #-check all -warn all -u
        FFLAGS   += -g -traceback -O0 -fpe0 #-check all -warn all -u
        CFLAGS   += -g -fp-trap=common #-Wcheck
        CXXFLAGS += -g -fp-trap=common #-Wcheck
      else
        F90FLAGS += -g -debug inline-debug-info -O2 -ip -align array64byte -qopt-report=5 -qopt-report-phase=vec
        FFLAGS   += -g -debug inline-debug-info -O2 -ip -align array64byte -qopt-report=5 -qopt-report-phase=vec
        CFLAGS   += -g -debug inline-debug-info -O2 -ip -qopt-report=5 -qopt-report-phase=vec
        CXXFLAGS += -g -debug inline-debug-info -O2 -ip -qopt-report=5 -qopt-report-phase=vec
      endif
      ifdef GPROF
        F90FLAGS += -p
        FFLAGS   += -p
        CFLAGS   += -p
        CXXFLAGS += -p
      endif
    endif

    ifeq ($(_comp),Intel16)
      ifndef NDEBUG
        F90FLAGS += -g -traceback -O0 -fpe0 #-check all -warn all -u 
        FFLAGS   += -g -traceback -O0 -fpe0 #-check all -warn all -u 
        CFLAGS   += -g -fp-trap=common #-Wcheck
        CXXFLAGS += -g -fp-trap=common #-Wcheck
      else
        F90FLAGS += -g -debug inline-debug-info -O2 -ip -align array64byte -qopt-report=5 -qopt-report-phase=vec
        FFLAGS   += -g -debug inline-debug-info -O2 -ip -align array64byte -qopt-report=5 -qopt-report-phase=vec
        CFLAGS   += -g -debug inline-debug-info -O2 -ip -qopt-report=5 -qopt-report-phase=vec
        CXXFLAGS += -g -debug inline-debug-info -O2 -ip -qopt-report=5 -qopt-report-phase=vec
      endif
      ifdef GPROF
        F90FLAGS += -pg
        FFLAGS   += -pg
        CFLAGS   += -pg
        CXXFLAGS += -pg
      endif
    endif

    ifeq ($(_comp),Intel15)
      ifndef NDEBUG
        F90FLAGS += -g -traceback -O0 -fpe0 #-check all -warn all -u 
        FFLAGS   += -g -traceback -O0 -fpe0 #-check all -warn all -u 
        CFLAGS   += -g -fp-trap=common #-Wcheck
        CXXFLAGS += -g -fp-trap=common #-Wcheck
      else
        F90FLAGS += -g -debug inline-debug-info -O2 -ip -align array64byte -qopt-report=5 -qopt-report-phase=vec
        FFLAGS   += -g -debug inline-debug-info -O2 -ip -align array64byte -qopt-report=5 -qopt-report-phase=vec
        CFLAGS   += -g -debug inline-debug-info -O2 -ip -qopt-report=5 -qopt-report-phase=vec
        CXXFLAGS += -g -debug inline-debug-info -O2 -ip -qopt-report=5 -qopt-report-phase=vec
      endif
      ifdef GPROF
        F90FLAGS += -pg
        FFLAGS   += -pg
        CFLAGS   += -pg
        CXXFLAGS += -pg
      endif
    endif

    ifeq ($(_comp),Intel14)
      ifndef NDEBUG
        F90FLAGS += -g -traceback -O0 #-check all -warn all -u 
        FFLAGS   += -g -traceback -O0 #-check all -warn all -u 
        CFLAGS   += -g #-Wcheck
        CXXFLAGS += -g #-Wcheck
      else
        F90FLAGS += -g -O2 -ip # -xHost # -fp-model source -vec-report6
        FFLAGS   += -g -O2 -ip # -xHost # -fp-model source 
        CFLAGS   += -g -O2 -ip # -xHost # -fp-model source 
        CXXFLAGS += -g -O2 -ip # -xHost # -fp-model source 
      endif
      ifdef GPROF
        F90FLAGS += -pg
        FFLAGS   += -pg
        CFLAGS   += -pg
        CXXFLAGS += -pg
      endif
    endif

    ifeq ($(_comp),Intel13)
      ifndef NDEBUG
        F90FLAGS += -g -traceback -O0 #-check all -warn all -u 
        FFLAGS   += -g -traceback -O0 #-check all -warn all -u 
        #CFLAGS   += -g #-Wcheck
      else
        ifdef INTEL_X86
	  F90FLAGS += -g -fast
	  FFLAGS   += -g -fast
	  CFLAGS   += -g -fast
	  CXXFLAGS += -g -fast
	else
          F90FLAGS += -g -O2 -ip -fp-model source #-xHost
          FFLAGS   += -g -O2 -ip -fp-model source #-xHost
          CFLAGS   += -g -O2 -ip -fp-model source #-xHost
          CXXFLAGS += -g -O2 -ip -fp-model source #-xHost
	endif
      endif
      ifdef GPROF
        F90FLAGS += -pg
      endif
#      F90FLAGS += -stand f95
#     FFLAGS += -stand f95
    endif
