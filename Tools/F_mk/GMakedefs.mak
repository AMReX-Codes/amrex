ARCH := $(shell uname)
UNAMEN := $(shell uname -n)
HOSTNAMEF := $(shell hostname -f)

ifndef HOST
  HOST := $(UNAMEN)
endif

FC       :=
F90      :=
CC       :=
CXX      :=
F90FLAGS :=
FFLAGS   :=
CFLAGS   :=
CXXFLAGS :=
FPP_DEFINES :=

FCOMP_VERSION :=

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
ifdef ACC
  acc_suffix 	:= .acc
endif
ifndef NDEBUG
  debug_suffix 	:= .debug
endif
ifdef TEST
  ifdef NDEBUG
    debug_suffix := .test
  endif
endif
ifdef MIC
  mic_suffix    := .mic
endif
ifdef SDC
  sdc_suffix 	:= .SDC
endif
ifdef ROSE
  rose_suffix   := .rose
endif
ifdef ZMQ
  zmq_suffix := .zmq
endif
ifdef HDF
  hdf_suffix := .hdf
endif

suf=$(ARCH).$(COMP)$(rose_suffix)$(debug_suffix)$(prof_suffix)$(mpi_suffix)$(omp_suffix)$(acc_suffix)$(mic_suffix)$(sdc_suffix)$(zmq_suffix)$(hdf_suffix)

sources     =
fsources    =
f90sources  =
sf90sources =
F90sources  =
csources    =
cxxsources  =
libraries   =
xtr_libraries =
hypre_libraries =
mpi_libraries =
mpi_include_dir =
mpi_lib_dir =

# test for endianness to distinguish older IBM Mac's from newer Intel ones
endian := $(shell perl -MConfig -e '$$i=index $$Config{byteorder},1;if ($$i == 0){print "little"};')

ifeq ($(ARCH),Darwin)
  ifeq ($(findstring little,$(endian)),little)
    CPPFLAGS += -DBL_$(ARCH)_little
  else
    CPPFLAGS += -DBL_$(ARCH)
  endif
else
  CPPFLAGS += -DBL_$(ARCH)
endif

CPPFLAGS += -DFORTRAN_BOXLIB
FPPFLAGS += -DFORTRAN_BOXLIB

ifdef TEST
  CPPFLAGS += -DBL_TESTING
endif

ifndef NDEBUG
  CPPFLAGS += -DDEBUG
endif

ifdef CUDA
  CPPFLAGS += -DCUDA
  FPPFLAGS += -DCUDA
endif

F_C_LINK := UNDERSCORE

odir=.
mdir=.
tdir=.

tdir = t/$(suf)
odir = $(tdir)/o
mdir = $(tdir)/m
hdir = t/html


ifeq ($(COMP),g95)
  include $(AMREX_HOME)/Tools/F_mk/comps/g95.mak
endif

ifeq ($(findstring gfortran, $(COMP)), gfortran)
  include $(AMREX_HOME)/Tools/F_mk/comps/gfortran.mak
endif

ifeq ($(COMP),xlf)
  include $(AMREX_HOME)/Tools/F_mk/comps/xlf.mak
endif

ifeq ($(ARCH),Darwin)
  ifeq ($(COMP),IBM)
    include $(AMREX_HOME)/Tools/F_mk/comps/Darwin_ibm.mak
  endif

  ifeq ($(COMP),Intel)
    include $(AMREX_HOME)/Tools/F_mk/comps/Darwin_intel.mak
  endif
endif

ifeq ($(ARCH),FreeBSD)
endif

ifeq ($(ARCH),Linux)
  ifeq ($(COMP),Cray)
    include $(AMREX_HOME)/Tools/F_mk/comps/Linux_cray.mak
  endif

  ifeq ($(COMP),PGI)
    include $(AMREX_HOME)/Tools/F_mk/comps/Linux_pgi.mak
  endif

  ifeq ($(COMP),SunStudio)
    include $(AMREX_HOME)/Tools/F_mk/comps/Linux_sunstudio.mak
  endif

  ifeq ($(COMP),PathScale)
    include $(AMREX_HOME)/Tools/F_mk/comps/Linux_pathscale.mak
  endif

  ifeq ($(COMP),Intel)
    include $(AMREX_HOME)/Tools/F_mk/comps/Linux_intel.mak
  endif

  # Gottingen machines
  ifeq ($(HOST),hicegate0)
    include $(AMREX_HOME)/Tools/F_mk/comps/Linux_intel.mak
  endif

  ifeq ($(COMP),NAG)
    include $(AMREX_HOME)/Tools/F_mk/comps/Linux_nag.mak
  endif

  ifeq ($(COMP),Lahey)
    include $(AMREX_HOME)/Tools/F_mk/comps/Linux_lahey.mak
  endif

  ifeq ($(findstring summit, $(HOST)), summit)
    override CPP_ARGS := -E
  endif
endif

ifeq ($(ARCH),AIX)
  include $(AMREX_HOME)/Tools/F_mk/comps/aix.mak
endif

ifeq ($(findstring mira, $(HOSTNAMEF)), mira)
  include $(AMREX_HOME)/Tools/F_mk/comps/bgq.mak
endif

ifeq ($(strip $(F90)),)
   $(error "COMP=$(COMP) is not supported")   
endif

ifdef ROSE
ifeq ($(strip $(ROSECOMP)),)
   $(error "ROSECOMP is not defined")   
endif
endif

ifdef MPI
  include $(AMREX_HOME)/Tools/F_mk/GMakeMPI.mak
endif

ifdef mpi_include_dir
  fpp_flags += -I$(mpi_include_dir)
endif

ifdef mpi_lib_dir
  fld_flags += -L$(mpi_lib_dir)
endif

f_includes = $(addprefix -I , $(FINCLUDE_LOCATIONS))
c_includes = $(addprefix -I , $(INCLUDE_LOCATIONS))

# MODDEP is for .f90, .f, and .F90.  
# MKDEP is for c
MODDEP  :=  $(AMREX_HOME)/Tools/F_scripts/dep.py
MKDEP   :=  $(AMREX_HOME)/Tools/F_scripts/mkdep.pl
F90DOC  :=  $(AMREX_HOME)/Tools/F_scripts/f90doc/f90doc

FPPFLAGS += $(fpp_flags) $(addprefix -I, $(FINCLUDE_LOCATIONS))
LDFLAGS  += $(fld_flags)
libraries += $(hypre_libraries) $(mpi_libraries) $(xtr_libraries)

CPPFLAGS += -DBL_FORT_USE_$(F_C_LINK) $(addprefix -I, $(INCLUDE_LOCATIONS))

objects = $(addprefix $(odir)/,       \
	$(sort $(f90sources:.f90=.o)) \
	$(sort $(sf90sources:.f90=.o)) \
	$(sort $(F90sources:.F90=.o)) \
	$(sort $(fsources:.f=.o))     \
	$(sort $(csources:.c=.o))     \
	$(sort $(cxxsources:.cpp=.o))     \
	)
sources =                     \
	$(sort $(f90sources)) \
        $(sort $(F90sources)) \
	$(sort $(fsources)  ) \
	$(sort $(csources)  ) \
	$(sort $(cxxsources)  )

html_sources = $(addprefix $(hdir)/,     \
	$(sort $(f90sources:.f90=.html)) \
	$(sort $(F90sources:.F90=.html)) \
	$(sort $(fsources:.f=.html))     \
	)

pnames = $(addsuffix .$(suf).exe, $(basename $(programs)))

ifeq ($(USE_CCACHE),TRUE)
  CCACHE = ccache
else
  CCACHE =
endif

ifeq ($(USE_F90CACHE),TRUE)
  F90CACHE = f90cache
else
  F90CACHE =
endif

ifndef ROSE
   COMPILE.f   = $(F90CACHE) $(FC)  $(FFLAGS)   $(FPPFLAGS) $(TARGET_ARCH) -c
   COMPILE.f90 = $(F90CACHE) $(F90) $(F90FLAGS) $(FPPFLAGS) $(TARGET_ARCH) -c
   COMPILE.c   = $(CCACHE)   $(CC)  $(CFLAGS)   $(CPPFLAGS) $(TARGET_ARCH) -c
   COMPILE.cc  = $(CCACHE)   $(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c
else
   COMPILE.cc  = $(ROSECOMP) $(CXXFLAGS) $(ROSEFLAGS) $(CPPFLAGS) -c
   COMPILE.c   = $(ROSECOMP) $(CFLAGS)   $(ROSEFLAGS) $(CPPFLAGS) -c
   COMPILE.f   = $(ROSECOMP) $(FFLAGS)   $(ROSEFLAGS) $(FPPFLAGS) -c
   COMPILE.f90 = $(ROSECOMP) $(F90FLAGS) $(ROSEFLAGS) $(FPPFLAGS) -c
endif

LINK.f      = $(FC)  $(FFLAGS) $(FPPFLAGS) $(LDFLAGS) $(TARGET_ARCH)
LINK.f90    = $(F90) $(F90FLAGS) $(FPPFLAGS) $(LDFLAGS) $(TARGET_ARCH)

# some pretty printing stuff
bold=`tput bold`
normal=`tput sgr0`

