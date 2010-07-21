ARCH := $(shell uname)
UNAMEN := $(shell uname -n)
HOSTNAMEF := $(shell hostname -f)

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

ifdef USE_HPCTOOLKIT
  hpc_suffix    := .hpc
endif
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

suf=$(ARCH).$(COMP)$(debug_suffix)$(prof_suffix)$(mpi_suffix)$(omp_suffix)$(hpc_suffix)

sources     =
fsources    =
f90sources  =
sf90sources  =
csources    =
libraries   =
xtr_libraries =
hypre_libraries =
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
  include '$(FPARALLEL)/mk/comps/g95.mak'
endif

# Note, we need a recent gfortran 4.2 build to compile --
# there are still runtime issues.
# to compile mt19937ar.f90, we need -fno-range-check, since
# that routine relies on overflows when doing initializations
ifeq ($(COMP),gfortran)
  include '$(FPARALLEL)/mk/comps/gfortran.mak'
endif

ifeq ($(COMP),xlf)
  include '$(FPARALLEL)/mk/comps/xlf.mak'
endif

ifeq ($(ARCH),Darwin)
  ifeq ($(COMP),IBM)
    include '$(FPARALLEL)/mk/comps/Darwin_ibm.mak'
  endif
endif

ifeq ($(ARCH),FreeBSD)
endif

ifeq ($(ARCH),Linux)
  ifeq ($(COMP),catamount)
    include '$(FPARALLEL)/mk/comps/Linux_catamount.mak'
  endif

  ifeq ($(COMP),xt4)
    include '$(FPARALLEL)/mk/comps/Linux_xt4.mak'
  endif

  ifeq ($(COMP),PGI)
    include '$(FPARALLEL)/mk/comps/Linux_pgi.mak'
  endif

  ifeq ($(COMP),SunStudio)
    include '$(FPARALLEL)/mk/comps/Linux_sunstudio.mak'
  endif

  ifeq ($(COMP),PathScale)
    include '$(FPARALLEL)/mk/comps/Linux_pathscale.mak'
  endif

  ifeq ($(COMP),Intel)
    include $(FPARALLEL)/mk/comps/Linux_intel.mak
  endif

  ifeq ($(COMP),NAG)
    include '$(FPARALLEL)/mk/comps/Linux_nag.mak'
  endif

  ifeq ($(COMP),Lahey)
    include '$(FPARALLEL)/mk/comps/Linux_lahey.mak'
  endif
endif

ifeq ($(ARCH),CRAYX1)
  include '$(FPARALLEL)/mk/comps/crayx1.mak'
endif

ifeq ($(ARCH),AIX)
  include '$(FPARALLEL)/mk/comps/aix.mak'
endif

ifeq ($(ARCH),IRIX64)
  include '$(FPARALLEL)/mk/comps/irix64.mak'
endif

ifeq ($(ARCH),OSF1)
  include '$(FPARALLEL)/mk/comps/osf1.mak'
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
libraries += $(hypre_libraries) $(mpi_libraries) $(xtr_libraries)

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
