ARCH := $(shell uname)

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
csources    =


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
  F90FLAGS += -std=f95
  FFLAGS   += -std=f95
  F90FLAGS += -fmod=$(mdir)
  FFLAGS   += -fmod=$(mdir)
  F90FLAGS += -Wall 
  FFLAGS += -Wall 
  CFLAGS += -Wall
  ifdef NDEBUG
    F90FLAGS += -O
    FFLAGS += -O
    CFLAGS += -O
  else
    F90FLAGS += -g
    F90FLAGS += -fbounds-check
    FFLAGS += -g
    FFLAGS += -fbounds-check
    CFLAGS += -g
  endif
  ifdef mpi_include
    fpp_flags += -I $(mpi_include)
  endif
  ifdef mpi_lib
    fld_flags += -L $(mpi_lib)
  endif
endif

# ALL GFORTRANS, which probably don't work, are the same
ifeq ($(COMP),gfortran)
  FC := gfortran
  F90 := gfortran
  CC := gcc
  F90FLAGS += -J$(mdir) -I $(mdir)
  FFLAGS   += -J$(mdir) -I $(mdir)
  CFLAGS += -Wall
  ifdef NDEBUG
    F90FLAGS += -O
    FFLAGS += -O
    CFLAGS += -O
  else
    F90FLAGS += -g 
    F90FLAGS += -fbounds-check
    FFLAGS += -g 
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
  ifdef MPI
    ifndef MPICHHOME
      $(error "MPI needed but MPIHOME not defined")
    endif
    mpi_home = $(MPICHHOME)
    mpi_include = $(mpi_home)/include
    mpi_lib = $(mpi_home)/lib
    mpi_libraries = -lmpich -lfmpich
  endif
  ifeq ($(COMP),PathScale)
    FC = pathf90
    F90 = pathf90
    FFLAGS += -nog77mangle
    F90FLAGS += -nog77mangle
    CC  = pathcc
    ifndef NDEBUG
      F90FLAGS += -g
      FFLAGS += -g
      CFLAGS += -g
#     F90FLAGS += -C
#     FFLAGS += -C
    else
      F90FLAGS += -Ofast
      FFLAGS += -Ofast
      CFLAGS += -Ofast
    endif
    LDFLAGS += -static
    CPPFLAGS += -DBL_HAS_SECOND
  endif
  ifeq ($(COMP),Intel)
    _unamem := $(shell uname -m)
    ifeq ($(_unamem),ia64)
      _ifc := efc
      _icc := icc
    else
      _ifc := ifc
      _icc := icc
    endif
    _ifc_version := $(shell $(_ifc) -V 2>&1 | grep 'Version')
    _icc_version := $(shell $(_icc) -V 2>&1 | grep 'Version')
    ifeq ($(findstring Version 8, $(_ifc_version)), Version 8)
      _ifc  := ifort
      _comp := Intel8
    else
    ifeq ($(findstring Version 7.1, $(_ifc_version)), Version 7.1)
      $(error "VERSION 7.1 of IFC Will Not Work")
    else
    ifeq ($(findstring Version 7.0, $(_ifc_version)), Version 7.0)
      _comp := Intel7
    else
      $(errorr "$(_ifc_version) of IFC will not work")
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
    FFLAGS   += -I $(mdir)
    F90FLAGS += -I $(mdir)
    ifdef OMP
      FFLAGS   += -openmp -fpp2
      F90FLAGS += -openmp -fpp2
    endif
    ifeq ($(_comp),Intel8)
      ifndef NDEBUG
        F90FLAGS += -g -traceback
        FFLAGS   += -g -traceback
        F90FLAGS += -check all
        FFLAGS   += -check all
        F90FLAGS += -check noshape -check nopointer
        FFLAGS   += -check noshape -check nopointer
        CFLAGS   += -g -Wcheck
	LDFLAGS  += -Bstatic
      else
        F90FLAGS += -O3
        FFLAGS += -O3
        CFLAGS += -O3
#       F90FLAGS += -ip
#       FFLAGS += -ip
#       CFLAGS += -ip
#       F90FLAGS += -ipo
#       FFLAGS += -ipo
#       CFLAGS += -ipo
#       LDFLAGS += -static
      endif
      ifdef PROF
        F90FLAGS += -pg
      endif
      # F90FLAGS += -stand f95
      # FFLAGS += -stand f95
    endif
    ifeq ($(_comp),Intel7)
      ifdef NDEBUG
        ifdef OMP
          CFLAGS += -O0
          FFLAGS += -O0
          F90FLAGS += -O0
	  F90FLAGS += -parallel -par_report3 -openmp_report2
        else
          CFLAGS += -O3
#         FFLAGS += -g
          FFLAGS += -O3
#         F90FLAGS += -g
          F90FLAGS += -O3
          ifndef PROF
#           CFLAGS += -ip
            CFLAGS += -ipo
#           FFLAGS += -ip
            FFLAGS += -ipo
#	    F90FLAGS += -ip
            F90FLAGS += -ipo
          endif
        endif
      else
#       F90FLAGS += -CU -CV -CS -CA -CB
        ifdef MPI
          FFLAGS   +=
          F90FLAGS +=
          CFLAGS   += 
        else
          CFLAGS += -g
          FFLAGS += -g
          F90FLAGS += -g
          ifdef OMP
            FFLAGS   +=     -CV -CS -CA -CB
            F90FLAGS +=     -CV -CS     -CB
          else
            FFLAGS   += -CB -CU
#           FFLAGS   +=         -CV -CS -CA
#           FFLAGS   +=                 -CA
            F90FLAGS += -CB -CU -CV -CS
#           F90FLAGS +=                 -CA
          endif
        endif
      endif

      ifdef PROF
        F90FLAGS += -pg
      endif

      ifdef MPI
        mpi_libraries += -lPEPCF90
      endif
      fld_flags  += -Vaxlib
    endif
    ifdef mpi_include
      fpp_flags += -I $(mpi_include)
    endif
    ifdef mpi_lib
      fld_flags += -L $(mpi_lib)
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
    ifdef PROF
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
    CFLAGS += -Wall
    ifdef NDEBUG
      FFLAGS += --tpp --prefetch 2 --nap --nchk -O --npca --nsav --ntrace
      F90FLAGS += --tpp --prefetch 2 --nap --nchk -O --npca --nsav --ntrace
    else
      FFLAGS   += -g --pca --nsav       --ap # --chk aesu # --chkglobal
      F90FLAGS += -g --pca --nsav --f95 --ap --chk aes  # --chkglobal
    endif
    ifdef mpi_include
      fpp_flags += -I $(mpi_include)
    endif
    ifdef mpi_lib
      fld_flags += -L $(mpi_lib)
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
  ifdef MPI
    F90 = mpxlf95$(rsuf)
    FC  = mpxlf$(rsuf)
  else
    F90 = xlf95$(rsuf)
    FC  = xlf$(rsuf)
  endif
  F90FLAGS += -qsuffix=f=f90 -qnosave
  F90FLAGS += -qmoddir=$(mdir)
  F90FLAGS += -I$(mdir)
  FC  += -qnosave
  ifdef NDEBUG
    ifdef PROF
      FFLAGS += -O2
      F90FLAGS += -O2
    else
      #FFLAGS += -O3 -qstrict -qtune=pwr3 -qarch=pwr3 -qipa -qhot
      FFLAGS += -O4
      #FFLAGS += -O2 -qstrict
      #F90FLAGS += -O3 -qstrict -qtune=pwr3 -qarch=pwr3 -qipa -qhot
      #F90FLAGS += -O3 -qipa -qhot
      F90FLAGS += -O4
      #F90FLAGS += -O2 -qstrict
    endif
  else
#   FFLAGS += -O0
    FFLAGS += -g
    FFLAGS += -C
    FFLAGS += -qinitauto=FF
#   F90FLAGS += -O0
    F90FLAGS += -g
    F90FLAGS += -C
    F90FLAGS += -qlanglvl=95std
    F90FLAGS += -qinitauto=FF
  endif
  ifdef OMP
    FFLAGS += -qsmp=omp
    F90FLAGS += -qsmp=omp
  endif
  ifdef PROF
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
  ifdef MPI
    mpi_libraries = -lmpi
  endif
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

c_includes = $(addprefix --I , $(INCLUDE_LOCATIONS))

TCSORT  :=  $(FPARALLEL)/scripts/tcsort.pl
MODDEP  :=  $(FPARALLEL)/scripts/moddep.pl
MKDEP   :=  $(FPARALLEL)/scripts/mkdep.pl
F90DOC  :=  $(FPARALLEL)/scripts/f90doc/f90doc

FPPFLAGS += $(fpp_flags)
LDFLAGS  += $(fld_flags)
libraries += $(mpi_libraries)

CPPFLAGS += -DBL_FORT_USE_$(F_C_LINK) $(addprefix -I, $(INCLUDE_LOCATIONS))

objects = $(addprefix $(odir)/,       \
	$(sort $(f90sources:.f90=.o)) \
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
