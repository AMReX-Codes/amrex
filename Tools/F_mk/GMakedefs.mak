ARCH := $(shell uname)
UNAMEN := $(shell uname -n)

ifeq ($(ARCH),UNICOS/mp)
  ARCH := CRAYX1
endif

FC       :=
F90      :=
F90FLAGS :=
FFLAGS   :=
CFLAGS   :=

VPATH_LOCATIONS :=
INCLUDE_LOCATIONS :=

ifdef MPI
  mpi_suffix 	:= .mpi
endif
ifdef PROF
  prof_suffix 	:= .prof
endif
ifdef OMP
  omp_suffix 	:= .omp
endif
ifndef NDEBUG
  debug_suffix 	:= .debug
endif

suf=$(ARCH).$(COMP)$(debug_suffix)$(prof_suffix)$(mpi_suffix)$(omp_suffix)

sources     =
fsources    =
f90sources  =
sf90sources  =
csources    =
libraries   =
xtr_libraries =
mpi_libraries =
mpi_include_dir =
mpi_lib_dir =

CPPFLAGS += -DBL_$(ARCH)

F_C_LINK := UNDERSCORE

odir=.
mdir=.
tdir=.

tdir = t/$(suf)
odir = $(tdir)/o
mdir = $(tdir)/m
hdir = t/html

# ALL G95's are the same
ifeq ($(COMP),g95)
  F_C_LINK := DBL_UNDERSCORE
  FC := g95
  F90 := g95
  CC := gcc
# F90FLAGS += -std=f95 -fintrinsic-extensions=iargc,getarg
# FFLAGS   += -std=f95 -fintrinsic-extensions=iargc,getarg  
  F90FLAGS += -fmod=$(mdir) -I $(mdir)
  FFLAGS   += -fmod=$(mdir) -I $(mdir)
# F90FLAGS += -Wall 
# FFLAGS += -Wall 
# CFLAGS += -Wall
# FFLAGS += -ffloat-store
# F90FLAGS += -ffloat-store
  ifdef NDEBUG
    F90FLAGS += -O3 -ffloat-store
    FFLAGS += -O3 -ffloat-store
    CFLAGS += -O3
  else
    F90FLAGS += -g
    F90FLAGS += -fbounds-check
    F90FLAGS += -freal=nan
    FFLAGS += -g
    FFLAGS += -fbounds-check
    FFLAGS += -freal=nan
    CFLAGS += -g
  endif
  ifdef GPROF
    F90FLAGS += -pg
    FFLAGS += -pg
    CFLAGS += -pg
  endif
  ifeq ($(ARCH),Darwin)
    xtr_libraries += -lSystemStubs
  endif
endif

# Note, we need a recent gfortran 4.2 build to compile --
# there are still runtime issues.
# to compile mt19937ar.f90, we need -fno-range-check, since
# that routine relies on overflows when doing initializations
ifeq ($(COMP),gfortran)
  FC := gfortran
  F90 := gfortran
  CC := gcc
  F90FLAGS += -J$(mdir) -I $(mdir)
  FFLAGS   += -J$(mdir) -I $(mdir)
  CFLAGS += -Wall
  ifdef NDEBUG
    F90FLAGS += -O -fno-range-check
    FFLAGS += -O -fno-range-check
    CFLAGS += -O
  else
    F90FLAGS += -g -fno-range-check 
    F90FLAGS += -fbounds-check 
    FFLAGS += -g -fno-range-check 
    FFLAGS += -fbounds-check 
    CFLAGS += -g
  endif
endif

ifeq ($(ARCH),Darwin)
  ifeq ($(COMP),IBM)
    F_C_LINK := LOWERCASE
    FC := xlf
    F90 := xlf95
    CC  := xlc
    F90FLAGS += -qsuffix=f=f90 -qnosave
    F90FLAGS += -qmoddir=$(mdir)
    ifdef NDEBUG
      F90FLAGS += -O5 -qtune=auto -qarch=auto -qunroll=auto
      FFLAGS   += -O5 -qtune=auto -qarch=auto -qunroll=auto
      CFLAGS   += -O5 -Q=20 -qtune-auto -qarch=auto -qunroll=auto -qaltivec
    else
      F90FLAGS += -g -C
      FFLAGS   += -g -C
    endif
    F90FLAGS += -I$(mdir)
  endif
endif

ifeq ($(ARCH),FreeBSD)
endif

ifeq ($(ARCH),Linux)
  ifeq ($(COMP),catamount)
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
  endif

  ifeq ($(COMP),xt4)
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
      FFLAGS   += -g
      F90FLAGS += -g
    endif
  endif

  ifeq ($(COMP),PGI)
    CC  := pgcc
    FC  := pgf95
    F90 := pgf95
    FFLAGS   += -module $(mdir) -I$(mdir) 
    F90FLAGS += -module $(mdir) -I$(mdir)
    ifdef NDEBUG
      FFLAGS   += -O
      F90FLAGS += -O
      CFLAGS   += -O
    else
      FFLAGS   += -g
      F90FLAGS += -g
    endif
  endif

  ifeq ($(COMP),SunStudio)
    FC = f95
    F90 = f95
  endif

  ifeq ($(COMP),PathScale)
    FC = pathf95
    F90 = pathf95
    FFLAGS   += -module $(mdir) -I$(mdir) 
    F90FLAGS += -module $(mdir) -I$(mdir)
    CC  = pathcc

    ifeq ($(findstring atlas, $(UNAMEN)), atlas)
    endif

#   F_C_LINK := DBL_UNDERSCORE
    ifndef NDEBUG

      F90FLAGS += -g -fno-second-underscore
      FFLAGS += -g -fno-second-underscore
      CFLAGS += -g -fno-second-underscore
#     F90FLAGS += -C
#     FFLAGS += -C
    else
      F90FLAGS += -O -ipa #-fno-second-underscore
      FFLAGS   += -O -ipa #-fno-second-underscore
      CFLAGS   += -O -ipa #-fno-second-underscore
    endif
#   LDFLAGS += -static
    CPPFLAGS += -DBL_HAS_SECOND
  endif

  ifeq ($(COMP),Intel)
    _unamem := $(shell uname -m)
    _ifc := ifort
    _icc := icc 
    _ifc_version := $(shell $(_ifc) -V 2>&1 | grep 'Version')
    _icc_version := $(shell $(_icc) -V 2>&1 | grep 'Version')
    ifeq ($(findstring Version 10, $(_ifc_version)), Version 10)
      _ifc  := ifort
      _comp := Intel10
    else
    ifeq ($(findstring Version 9, $(_ifc_version)), Version 9)
      ifeq ($(findstring atlas, $(UNAMEN)), atlas)
        _ifc  := mpiifort
        _comp := Intel9
      else
        _ifc  := ifort
        _comp := Intel9
      endif
    else
    ifeq ($(findstring Version 8, $(_ifc_version)), Version 8)
      _ifc  := ifort
      _comp := Intel8
    else
      $(errorr "$(_ifc_version) of IFC is not supported")
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
#    F90FLAGS += -mp
    FFLAGS   += -I $(mdir)
    F90FLAGS += -I $(mdir)
    ifdef OMP
      FFLAGS   += -openmp -fpp2
      F90FLAGS += -openmp -fpp2
    endif
    ifeq ($(_comp),Intel10)
      ifndef NDEBUG
        F90FLAGS += -g -traceback -O0
        FFLAGS   += -g -traceback -O0
        F90FLAGS += -check all
        FFLAGS   += -check all
        #CFLAGS   += -g -Wcheck
      else
        ifdef INTEL_X86
	  F90FLAGS += -fast
	  FFLAGS += -fast
	  CFLAGS += -fast
	else
          F90FLAGS += -O3 -ip
          FFLAGS += -O3 -ip
          CFLAGS += -O3 -ip
	endif
      endif
      ifdef GPROF
        F90FLAGS += -pg
      endif
      F90FLAGS += -stand f95
#     FFLAGS += -stand f95
    endif
    ifeq ($(_comp),Intel9)
      ifndef NDEBUG
        F90FLAGS += -g -traceback -O0
        FFLAGS   += -g -traceback -O0
        F90FLAGS += -check all
        FFLAGS   += -check all
        #CFLAGS   += -g -Wcheck
      else
        ifdef INTEL_X86
	  F90FLAGS += -fast
	  FFLAGS += -fast
	  CFLAGS += -fast
	else
          F90FLAGS += -O3 -ip
          FFLAGS += -O3 -ip
          CFLAGS += -O3 -ip
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
      F90FLAGS += -stand f95
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
#	LDFLAGS  += -Bstatic
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
  endif
  ifeq ($(COMP),NAG)
    FC  = nf95
    F90 = nf95
    CC  = gcc
    F90FLAGS += -mdir $(mdir) -I $(mdir)
    FFLAGS   += -mdir $(mdir) -I $(mdir)
    FFLAGS   += -w=x77 -fixed
    CFLAGS += -Wall
#   F90FLAGS += -Oassumed=always_contig
    f2kcli_suf := _nag
    ifdef NDEBUG
      FFLAGS += -O4
      F90FLAGS += -O4
    else
      FFLAGS += -C=all
      F90FLAGS += -C=all
      FFLAGS += -g
      F90FLAGS += -g
      FFLAGS += -nan
      F90FLAGS += -nan
      FFLAGS += -gline
      F90FLAGS += -gline
    endif
    ifdef GPROF
      FFLAGS += -pg
      F90FLAGS += -pg
      CFLAGS += -pg
    endif
  endif
  ifeq ($(COMP),Lahey)
    FC  = lf95
    F90 = lf95
    CC  = gcc
    FFLAGS =
    F90FLAGS =
    F90FLAGS += -M $(mdir)
    FFLAGS += -M $(mdir)
    CFLAGS += -Wall
    ifdef NDEBUG
      FFLAGS += --tpp --prefetch 2 --nap --nchk -O --npca --nsav --ntrace
      F90FLAGS += --tpp --prefetch 2 --nap --nchk -O --npca --nsav --ntrace
    else
      FFLAGS   += -g --pca --nsav       --ap # --chk aesu # --chkglobal
      F90FLAGS += -g --pca --nsav --f95 --ap --chk aes  # --chkglobal
    endif
    ifdef OMP
      FFLAGS += --parallel --openmp 
      F90FLAGS += --parallel --openmp
    endif
  endif
endif

ifeq ($(ARCH),CRAYX1)
  COMP = cftn
  F90 := ftn
  FC  := ftn
  CC   := cc
  FFLAGS   =
  F90FLAGS =
  FFLAGS   += -p $(mdir)  -J $(mdir) -e m
  F90FLAGS += -p $(mdir)  -J $(mdir) -e m
  ifndef NDEBUG
    FFLAGS   += -g
    F90FLAGS += -g
  endif
  f2kcli_suf := _crayx1
endif

ifeq ($(ARCH),AIX)
  F_C_LINK := LOWERCASE
  COMP = xlf
  ifdef OMP
    rsuf := _r
  endif
  F90 = xlf95$(rsuf)
  FC  = xlf95$(rsuf)
  F90FLAGS += -qnosave -qmoddir=$(mdir) -I$(mdir) -qsuffix=f=f90
  FFLAGS   += -qnosave -qmoddir=$(mdir) -I$(mdir) -qsuffix=f=f -qfixed=72
  ifdef NDEBUG
    ifdef GPROF
      FFLAGS   += -O2
      F90FLAGS += -O2
    else
      FFLAGS   += -O3 -qstrict -qtune=auto -qarch=auto -qcache=auto -NS5000
      F90FLAGS += -O3 -qstrict -qtune=auto -qarch=auto -qcache=auto -NS5000
    endif
  else
    FFLAGS += -g
    FFLAGS += -C
    FFLAGS += -qinitauto=FF 
    FFLAGS += -qlanglvl=95std
    F90FLAGS += -g
    F90FLAGS += -C
    F90FLAGS += -qinitauto=FF
    F90FLAGS += -qlanglvl=95std
  endif
  #
  # You might need the following on seaborg:
  #
  # LDFLAGS += -bmaxdata:0x80000000

  ifdef OMP
    FFLAGS += -qsmp=omp
    F90FLAGS += -qsmp=omp
  endif
  ifdef GPROF
    FFLAGS += -pg
    F90FLAGS += -pg
  endif
endif

ifeq ($(ARCH),IRIX64)
  COMP = f90
  F90  = f90
  FC   = f90
  CC   = cc
  tdir = .
  odir = .
  mdir = .
  ifdef NDEBUG
#   FFLAGS += -O
#   F90FLAGS += -O
      FFLAGS += -Ofast
      F90FLAGS += -Ofast
      LDFLAGS += -Ofast
  else
#   FFLAGS += -C
    FFLAGS += -DEBUG:verbose_runtime=ON
    FFLAGS += -DEBUG:subscript_check=ON
#   FFLAGS += -DEBUG:trap_uninitialized=ON
    FFLAGS += -g
    FFLAGS += -ansi
#   F90FLAGS += -C
    F90FLAGS += -DEBUG:verbose_runtime=ON
    F90FLAGS += -DEBUG:conform_check=ON
    F90FLAGS += -DEBUG:subscript_check=ON
#   F90FLAGS += -DEBUG:trap_uninitialized=ON
    F90FLAGS += -g
    F90FLAGS += -ansi
  endif
  F90FLAGS += -64
  FFLAGS += -64
  CFLAGS += -64
  ifdef OMP
    F90FLAGS += -mp
    FFLAGS += -mp
  endif
endif

ifeq ($(ARCH),OSF1)
  COMP = f90
  F90 = f90
  FC  = f90
  ifdef DEBUG
    FFLAGS += -g   -check bounds
    F90FLAGS += -g  -check bounds
  else
    FFLAGS += -fast -inline all
    F90FLAGS += -fast -inline all
  endif
  ifdef OMP
    FFLAGS += -omp
    F90FLAGS += -omp
    ifdef DEBUG
      FFLAGS += -check omp_bindings
      F90FLAGS += -check omp_bindings
    endif
  endif
endif

ifdef MPI
  include $(FPARALLEL)/mk/GMakeMPI.mak
endif

ifdef mpi_include_dir
  fpp_flags += -I $(mpi_include_dir)
endif

ifdef mpi_lib_dir
  fld_flags += -L $(mpi_lib_dir)
endif

f_includes = $(addprefix -I , $(FINCLUDE_LOCATIONS))
c_includes = $(addprefix -I , $(INCLUDE_LOCATIONS))

TCSORT  :=  $(FPARALLEL)/scripts/tcsort.pl
MODDEP  :=  $(FPARALLEL)/scripts/moddep.pl
MKDEP   :=  $(FPARALLEL)/scripts/mkdep.pl
F90DOC  :=  $(FPARALLEL)/scripts/f90doc/f90doc

FPPFLAGS += $(fpp_flags) $(f_includes)
LDFLAGS  += $(fld_flags)
libraries += $(mpi_libraries) $(xtr_libraries)

CPPFLAGS += -DBL_FORT_USE_$(F_C_LINK) $(addprefix -I, $(INCLUDE_LOCATIONS))

objects = $(addprefix $(odir)/,       \
	$(sort $(f90sources:.f90=.o)) \
	$(sort $(sf90sources:.f90=.o)) \
	$(sort $(fsources:.f=.o))     \
	$(sort $(csources:.c=.o))     \
	)
sources =                     \
	$(sort $(f90sources)) \
	$(sort $(fsources)  ) \
	$(sort $(csources)  )

html_sources = $(addprefix $(hdir)/,     \
	$(sort $(f90sources:.f90=.html)) \
	$(sort $(fsources:.f=.html))     \
	)

pnames = $(addsuffix .$(suf).exe, $(basename $(programs)))

COMPILE.f   = $(FC)  $(FFLAGS) $(FPPFLAGS) $(TARGET_ARCH) -c
COMPILE.f90 = $(F90) $(F90FLAGS) $(FPPFLAGS) $(TARGET_ARCH) -c

LINK.f      = $(FC)  $(FFLAGS) $(FPPFLAGS) $(LDFLAGS) $(TARGET_ARCH)
LINK.f90    = $(F90) $(F90FLAGS) $(FPPFLAGS) $(LDFLAGS) $(TARGET_ARCH)


