    _ifc_version := $(shell ifort -V 2>&1 | grep 'Version')
    ifeq ($(findstring Version 15, $(_ifc_version)), Version 15)
        _comp := Intel15
    else ifeq ($(findstring Version 14, $(_ifc_version)), Version 14)
        _comp := Intel14
    else ifeq ($(findstring Version 13, $(_ifc_version)), Version 13)
        _comp := Intel13
    else ifeq ($(findstring Version 12, $(_ifc_version)), Version 12)
        _comp := Intel12
    else ifeq ($(findstring Version 11, $(_ifc_version)), Version 11)
        _comp := Intel11
    else ifeq ($(findstring Version 10, $(_ifc_version)), Version 10)
        _comp := Intel10
    else
      $(error "$(_ifc_version) of IFC is not supported")
    endif

    F90 := ifort
    FC  := ifort
    CC  := icc
    CXX := icpc

    FFLAGS   += -module $(mdir) -I $(mdir)
    F90FLAGS += -module $(mdir) -I $(mdir) -cxxlib
    CFLAGS   += -std=c99

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
      FFLAGS   += -openmp
      F90FLAGS += -openmp
      CFLAGS   += -openmp
      CXXFLAGS += -openmp
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
        F90FLAGS += -O2 -ip # -xHost # -fp-model source -vec-report6
        FFLAGS   += -O2 -ip # -xHost # -fp-model source 
        CFLAGS   += -O2 -ip # -xHost # -fp-model source 
        CXXFLAGS += -O2 -ip # -xHost # -fp-model source 
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
	  F90FLAGS += -fast
	  FFLAGS += -fast
	  CFLAGS += -fast
	  CXXFLAGS += -fast
	else
          F90FLAGS += -O2 -ip -fp-model source #-xHost
          FFLAGS   += -O2 -ip -fp-model source #-xHost
          CFLAGS   += -O2 -ip -fp-model source #-xHost
          CXXFLAGS += -O2 -ip -fp-model source #-xHost
	endif
      endif
      ifdef GPROF
        F90FLAGS += -pg
      endif
#      F90FLAGS += -stand f95
#     FFLAGS += -stand f95
    endif

    ifeq ($(_comp),Intel12)
      ifndef NDEBUG
        F90FLAGS += -g -traceback -O0 #-check all -warn all -u 
        FFLAGS   += -g -traceback -O0 #-check all -warn all -u 
        #CFLAGS   += -g -Wcheck
      else
        ifdef INTEL_X86
	  F90FLAGS += -fast
	  FFLAGS += -fast
	  CFLAGS += -fast
	  CXXFLAGS += -fast
	else
          F90FLAGS += -O2 -ip -fp-model source #-xHost
          FFLAGS   += -O2 -ip -fp-model source #-xHost
          CFLAGS   += -O2 -ip -fp-model source #-xHost
          CXXFLAGS += -O2 -ip -fp-model source #-xHost
	endif
      endif
      ifdef GPROF
        F90FLAGS += -pg
      endif
#      F90FLAGS += -stand f95
#     FFLAGS += -stand f95
    endif

    ifeq ($(_comp),Intel11)
      ifndef NDEBUG
        F90FLAGS += -g -traceback -O0 -check all -warn all -u
        FFLAGS   += -g -traceback -O0 -check all -warn all -u
        #CFLAGS   += -g -Wcheck
      else
        ifdef INTEL_X86
          F90FLAGS += -fast
          FFLAGS += -fast
          CFLAGS += -fast
          CXXFLAGS += -fast
        else
          F90FLAGS += -O3 -ip -mp1# -fltconsistency
          FFLAGS += -O3 -ip -mp1# -fltconsistency
          CFLAGS += -O3 -ip -mp1
          CXXFLAGS += -O3 -ip -mp1
        endif
      endif
      ifdef GPROF
        F90FLAGS += -pg
      endif
#      F90FLAGS += -stand f95
#     FFLAGS += -stand f95
    endif

    ifeq ($(_comp),Intel10)
      ifndef NDEBUG
        F90FLAGS += -g -traceback -O0
        FFLAGS   += -g -traceback -O0
        F90FLAGS += -check all -warn all -u 
	#F90FLAGS += -ftrapuv
        FFLAGS   += -check all -warn all -u 
	#FFLAGS += -ftrapuv
        #CFLAGS   += -g -Wcheck
      else
        ifdef INTEL_X86
	  F90FLAGS += -fast
	  FFLAGS += -fast
	  CFLAGS += -fast
	  CXXFLAGS += -fast
	else
         # A. Donev added this to make moderately optimized executables:         
         ifndef BL_FAST_COMP
           F90FLAGS += -O3 -ip -mp1 -fltconsistency 
           FFLAGS += -O3 -ip -mp1 -fltconsistency
           CFLAGS += -O3 -ip -mp1 -fltconsistency
           CXXFLAGS += -O3 -ip -mp1 -fltconsistency
         else
           # Fast compiles and fast-enough runs:
           F90FLAGS += -O2 -mp1 -fltconsistency
           FFLAGS += -O2 -mp1 -fltconsistency
           CFLAGS += -O2 -mp1
           CXXFLAGS += -O2 -mp1
         endif        
	endif
      endif
      ifdef GPROF
        F90FLAGS += -pg
      endif
#      F90FLAGS += -stand f95
#     FFLAGS += -stand f95
    endif

