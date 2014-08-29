    _unamem := $(shell uname -m)
    _ifc := ifort
    _icc := icc
    _ifc_version := $(shell $(_ifc) -V 2>&1 | grep 'Version')
    _icc_version := $(shell $(_icc) -V 2>&1 | grep 'Version')
    ifeq ($(findstring Version 12, $(_ifc_version)), Version 12)
        _ifc  := ifort
        _comp := Intel12
    else
    ifeq ($(findstring Version 11, $(_ifc_version)), Version 11)
        _ifc  := ifort
        _comp := Intel11
    else
    ifeq ($(findstring Version 10, $(_ifc_version)), Version 10)
        _ifc  := ifort
        _comp := Intel10
    else
    ifeq ($(findstring Version 9, $(_ifc_version)), Version 9)
        _ifc  := ifort
        _comp := Intel9
    else
    ifeq ($(findstring Version 8, $(_ifc_version)), Version 8)
      _ifc  := ifort
      _comp := Intel8
    else
      $(error "$(_ifc_version) of IFC is not supported")
    endif
    endif
    endif
    endif
    endif
#   _ifc += -auto
    F90 := $(_ifc)
    FC  := $(_ifc)
    CC  := $(_icc)
    FFLAGS   =
    F90FLAGS =
    CFLAGS   =
    FFLAGS   += -module $(mdir)
    F90FLAGS += -module $(mdir)
    F90FLAGS += -mp1 -fltconsistency
    FFLAGS   += -I $(mdir)
    F90FLAGS += -I $(mdir)
    CFLAGS   += -std=c99
    ifdef OMP
      FFLAGS   += -openmp -openmp-report2
      F90FLAGS += -openmp -openmp-report2
    endif
    ifeq ($(_comp),Intel12)
      ifndef NDEBUG
        F90FLAGS += -g -traceback -O0 -check all -warn all -u
        FFLAGS   += -g -traceback -O0 -check all -warn all -u
        #CFLAGS   += -g -Wcheck
      else
        ifdef INTEL_X86
         F90FLAGS += -fast
         FFLAGS += -fast
         CFLAGS += -fast
       else
          F90FLAGS += -O3 -ip -mp1# -fltconsistency
          FFLAGS += -O3 -ip -mp1# -fltconsistency
          CFLAGS += -O3 -ip -mp1
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
        else
          F90FLAGS += -O3 -ip -mp1# -fltconsistency
          FFLAGS += -O3 -ip -mp1# -fltconsistency
          CFLAGS += -O3 -ip -mp1
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
       else
         # A. Donev added this to make moderately optimized executables:
         ifndef BL_FAST_COMP
           F90FLAGS += -O3 -ip -mp1 -fltconsistency
           FFLAGS += -O3 -ip -mp1 -fltconsistency
           CFLAGS += -O3 -ip -mp1 -fltconsistency
         else
           # Fast compiles and fast-enough runs:
           F90FLAGS += -O2 -mp1 -fltconsistency
           FFLAGS += -O2 -mp1 -fltconsistency
           CFLAGS += -O2 -mp1
         endif
       endif
      endif
      ifdef GPROF
        F90FLAGS += -pg
      endif
#      F90FLAGS += -stand f95
#     FFLAGS += -stand f95
    endif
    ifeq ($(_comp),Intel9)
      ifndef NDEBUG
        F90FLAGS += -g -traceback -O0
        FFLAGS   += -g -traceback -O0
        F90FLAGS += -check all -warn all -u
        FFLAGS   += -check all -warn all -u
        #CFLAGS   += -g -Wcheck
      else
        ifdef INTEL_X86
         F90FLAGS += -fast
         FFLAGS += -fast
         CFLAGS += -fast
       else
          F90FLAGS += -O -mp1
          FFLAGS += -O -mp1
          CFLAGS += -O -mp1
#  ifndef GPROF
          F90FLAGS += #-ipo
          FFLAGS += #-ipo
          CFLAGS += #-ipo
#  endif
       endif
      endif
      ifdef GPROF
        F90FLAGS += -pg
      endif
#      F90FLAGS += -stand f95
#     FFLAGS += -stand f95
    endif
    ifeq ($(_comp),Intel8)
      ifndef NDEBUG
        F90FLAGS += -g -traceback
        FFLAGS   += -g -traceback
        F90FLAGS += -check all
        FFLAGS   += -check all
#       F90FLAGS += -check noshape -check nopointer
#       FFLAGS   += -check noshape -check nopointer
        CFLAGS   += -g -Wcheck
#      LDFLAGS  += -Bstatic
      else
        ifdef INTEL_X86
         F90FLAGS += -fast
         FFLAGS += -fast
         CFLAGS += -fast
       else
          F90FLAGS += -O3
          FFLAGS += -O3
          CFLAGS += -O3
          ifndef GPROF
            F90FLAGS += #-ipo
            FFLAGS += #-ipo
            CFLAGS += #-ipo
          endif
       endif
#       LDFLAGS += -static
      endif
      ifdef GPROF
        F90FLAGS += -pg
      endif
      # F90FLAGS += -stand f95
      # FFLAGS += -stand f95
    endif
