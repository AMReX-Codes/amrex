NDEBUG :=
NDEBUG :=
MPI    :=
OMP    :=

MULTIGRID := t
#VARDEN    := t
EIGEN     := t

ARCH := $(shell uname)
COMP := PathScale
COMP := Lahey8
COMP := Intel7
COMP := Intel


FC       :=
F90      :=
F90FLAGS :=
FFLAGS   :=
CFLAGS   :=

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
csources    =
mpi_sources =
omp_sources =

ifndef MPI
  mpi_sources += parallel_stubs.f90
else
  mpi_sources += parallel.f90
  mpi_sources += mpi.f
endif

sources += main.f90
sources += t_main.f90
sources += t_stencil.f90

sources += interp.f90

ifdef VARDEN
sources += varden.f90
sources += macproject.f90
sources += hgproject.f90
sources += advance.f90
sources += estdt.f90
sources += initdata.f90
sources += mkflux.f90
sources += slope.f90
sources += bc.f90
sources += mkutrans.f90
sources += cvmg.f90
sources += setvelbc.f90
endif

ifdef EIGEN
sources += cc_eigen.f90
sources += t_eigen.f90
endif

ifdef MULTIGRID
sources += wrapper.f90
sources += cc_single.f90 cc_multi.f90
sources += nodal_multi.f90
sources += ml.f90
sources += mg.f90
sources += mg_smoother.f90
sources += itsol.f90
sources += sparse_solve.f90
sources += stencil.f90
sources += stencil_nodal.f90
sources += mlmg.f90
sources += bndry_reg.f90

sources += mg_cpp.f90
endif

sources += BoxLib.f90
sources += omp.f90
sources += f2kcli$(f2kcli_suf).f90
sources += mt19937ar.f90
sources += bl_space.f90  
sources += bl_error.f90  
sources += bl_types.f90  
sources += bl_utility.f90
sources += bl_constants.f90
sources += box.f90  
sources += knapsack.f90 
sources += layout.f90  
sources += boxarray.f90  
sources += mboxarray.f90
sources += box_util.f90
sources += fab.f90  
sources += multifab.f90
sources += fabio.f90
sources += plotfile.f90
sources += cluster.f90

sources += list_box.f90
sources += sort_box.f90
sources += vector_i.f90 
sources += sort_d.f90

sources += ppm_util.f90

#sources += interp.f90

#sources += sort_i.f90
ifdef MPI
sources += pingpong.f90
endif

# SPARSKIT files
ifdef MULTIGRID
sources += ilut.f iters.f sk_sup.f
endif

ifdef EIGEN
# LAPAC files
sources += dsygv.f
endif

sources += $(mpi_sources)
sources += $(omp_sources)

csources += fabio_c.c
csources += timer_c.c
csources += ppm_util_c.c

CPPFLAGS += -DBL_$(ARCH)
ifeq ($(ARCH),AIX)
CPPFLAGS += -DBL_FORT_USE_LOWERCASE
else
ifeq ($(COMP),g95)
CPPFLAGS += -DBL_FORT_USE_DBL_UNDERSCORE
else
CPPFLAGS += -DBL_FORT_USE_UNDERSCORE
endif
endif

odir=.
mdir=.
tdir=.
# hdir=./doc/html

tdir = t/$(suf)
odir = $(tdir)/o
mdir = $(tdir)/m

ifeq ($(ARCH),Linux)
  ifdef MPI
    mpi_home = /usr/local/mpi
    mpi_include = $(mpi_home)/include
    mpi_lib = $(mpi_home)/lib
    mpi_libraries = -lmpich
  endif
  ifeq ($(COMP),g95)
    FC := g95
    F90 := g95
    CC := gcc
    F90FLAGS += -fmod=$(mdir) -I $(mdir)
    FFLAGS   += -fmod=$(mdir) -I $(mdir)
    ifndef OMP
      omp_sources += omp_stubs.f90
    endif
    ifdef NDEBUG
      F90FLAGS += -O
      FFLAGS += -O
      CFLAGS += -O
    else
      F90FLAGS += -g -Wsurprising -fbounds-check
      F90FLAGS += -fbounds-check
      FFLAGS += -g -Wsurprising -fbounds-check
      FFLAGS += -fbounds-check
      CFLAGS += -g
    endif
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
    ifndef OMP
      omp_sources += omp_stubs.f90
    endif
    LDFLAGS += -static
    CPPFLAGS += -DBL_HAS_SECOND
  endif
  ifeq ($(COMP),Intel8)
    FC = ifort
    F90 = ifort
    CC = icc
    FFLAGS =
    F90FLAGS =
    ifndef NDEBUG
                                # Someday, can eliminate these	
      F90FLAGS += -g -check all -check nopointer -check noshape
      FFLAGS   += -g -check all -check nopointer -check noshape
      CFLAGS   += -g
    else
      F90FLAGS += -O3
      FFLAGS += -O3
      CFLAGS += -O3
#     F90FLAGS += -fast
#     FFLAGS += -fast
#     CFLAGS += -fast
    endif
    ifndef OMP
      omp_sources += omp_stubs.f90
    else
      F90FLAGS += -openmp
      FFLAGS   += -openmp
    endif
  endif
  ifeq ($(COMP),Intel)
#   ifdef MPI
#     F90=mpif90
#     FC= mpif77 -cm -w90
#   else
      F90 = ifc
      FC  = ifc
#     FC  = ifc -cm -w90
#   endif
    CC  = icc
    FFLAGS   =
    F90FLAGS =
#   FFLAGS += -auto
#   F90FLAGS += -auto
    F90FLAGS += -module $(mdir)
    F90FLAGS += -I $(mdir)
    ifdef NDEBUG
      ifdef OMP
        CFLAGS += -O0
        FFLAGS += -O0
        F90FLAGS += -O0
	F90FLAGS += -parallel -par_report3 -openmp_report2
      else
        CFLAGS += -O3
#       FFLAGS += -g
        FFLAGS += -O3
#       F90FLAGS += -g
        F90FLAGS += -O3
        ifndef PROF
#         CFLAGS += -ip
          CFLAGS += -ipo
#         FFLAGS += -ip
          FFLAGS += -ipo
#	  F90FLAGS += -ip
          F90FLAGS += -ipo
        endif
      endif
    else
#     F90FLAGS += -CU -CV -CS -CA -CB
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
          FFLAGS   += -CU -CV -CS -CA -CB
#         FFLAGS   +=             -CA
          F90FLAGS += -CU -CV -CS     -CB
#         F90FLAGS +=             -CA
        endif
      endif
    endif

    ifdef PROF
      F90FLAGS += -pg
    endif

    ifdef OMP
      FFLAGS   += -openmp -fpp2
      F90FLAGS += -openmp -fpp2
    else
      FFLAGS   += -openmp_stubs
      F90FLAGS += -openmp_stubs
    endif
    ifdef mpi_include
      fpp_flags += -I $(mpi_include)
    endif
    ifdef mpi_lib
      fld_flags += -L $(mpi_lib)
    endif
    ifdef MPI
      mpi_libraries += -lPEPCF90
    endif
    fld_flags  += -Vaxlib
  endif
  ifeq ($(COMP),NAG)
    FC  = nf95
    F90 = nf95
    CC  = gcc
    F90FLAGS += -mdir $(mdir) -I $(mdir)
    FFLAGS   += -mdir $(mdir) -I $(mdir)
    FFLAGS   += -w=x77 -fixed
#   F90FLAGS += -Oassumed=always_contig
    omp_sources += omp_stubs.f90
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
    ifdef NDEBUG
      FFLAGS += --tpp --prefetch 2 --nap --nchk -O --npca --nsav --ntrace
      F90FLAGS += --tpp --prefetch 2 --nap --nchk -O --npca --nsav --ntrace
    else
      FFLAGS   += -g --pca --nsav       --ap --chk aesu # --chkglobal
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
    else
      omp_sources += omp_stubs.f90
    endif
  endif
endif

ifeq ($(ARCH),AIX)
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
    FFLAGS += -DEBUG:trap_uninitialized=ON
    FFLAGS += -g
    FFLAGS += -ansi
#   F90FLAGS += -C
    F90FLAGS += -DEBUG:verbose_runtime=ON
    F90FLAGS += -DEBUG:conform_check=ON
    F90FLAGS += -DEBUG:subscript_check=ON
    F90FLAGS += -DEBUG:trap_uninitialized=ON
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
  ifndef OMP
    omp_sources += omp_stubs.f90
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

vpath %.f . ./SPARSKIT

TCSORT  :=  tcsort.pl
MODDEP  :=  moddep.pl
F90DOC  :=  f90doc


FPPFLAGS += $(fpp_flags)
LDFLAGS  += $(fld_flags)
libraries += $(mpi_libraries)

objects := $(addprefix $(odir)/, $(sources:.f90=.o) $(csources:.c=.o))
objects := $(objects:.f=.o)

# hsources := $(addprefix $(hdir)/, $(sources:.f90=.html))
# hsources := $(hsources:.f=.html)

main.$(suf).exe: $(objects)
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

# doc:	$(hsources)

COMPILE.f   = $(FC)  $(FFLAGS) $(FPPFLAGS) $(TARGET_ARCH) -c
COMPILE.f90 = $(F90) $(F90FLAGS) $(FPPFLAGS) $(TARGET_ARCH) -c

LINK.f      = $(FC)  $(FFLAGS) $(FPPFLAGS) $(LDFLAGS) $(TARGET_ARCH)
LINK.f90    = $(F90) $(F90FLAGS) $(FPPFLAGS) $(LDFLAGS) $(TARGET_ARCH)

${odir}/%.o: %.f
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	@if [ ! -d $(mdir) ]; then mkdir -p $(mdir); fi
	$(COMPILE.f) $(OUTPUT_OPTION) $<

${odir}/%.o: %.f90
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	@if [ ! -d $(mdir) ]; then mkdir -p $(mdir); fi
	$(COMPILE.f90) $(OUTPUT_OPTION) $<

${odir}/%.o: %.c
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	$(COMPILE.c) $(OUTPUT_OPTION) $<

ifeq ($(COMP),Intel)
$(odir)/dsygv.o: dsygv.f
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	@if [ ! -d $(mdir) ]; then mkdir -p $(mdir); fi
ifdef NDEBUG
ifdef PROF
	$(FC) -mp -O3      -c -o $(odir)/dsygv.o dsygv.f
else
	$(FC) -mp -O3 -ipo -c -o $(odir)/dsygv.o dsygv.f
endif
else
	$(FC) -g -mp -CB -CU -c -o $(odir)/dsygv.o dsygv.f
endif
endif
ifeq ($(COMP),g95)
	$(FC)              -c -o $(odir)/dsygv.o dsygv.f
endif

ifeq ($(COMP),Lahey)
$(odir)/dsygv.o: dsygv.f
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	@if [ ! -d $(mdir) ]; then mkdir -p $(mdir); fi
ifdef NDEBUG
	$(FC) -O -c -o $(odir)/dsygv.o dsygv.f
else
	$(FC) -g -c -o $(odir)/dsygv.o dsygv.f
endif
endif

ifdef MPI
$(odir)/mpi.o: mpi.f
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	@if [ ! -d $(mdir) ]; then mkdir -p $(mdir); fi
	$(FC) -c $(FPPFLAGS) -o $(odir)/mpi.o mpi.f
endif

# ${hdir}/%.html: %.f90
# 	@if [ ! -d $(hdir) ]; then mkdir -p $(hdir); fi
# 	$(F90DOC) $(OUTPUT_OPTION) $<

# ${hdir}/%.html: %.f
# 	@if [ ! -d $(hdir) ]; then mkdir -p $(hdir); fi
# 	$(F90DOC) $(OUTPUT_OPTION) $<

all.f90: $(sources)
	perl moddep.pl --tsort $(sources) | perl tcsort.pl > tmp.files
	cat `cat tmp.files` > all.f90

clean:
	$(RM) ./*.o ./*.mod $(mdir)/*.mod $(odir)/*.o *.exe *~
	$(RM) $(odir)/*.il
	$(RM) $(tdir)/f90.depends
	$(RM) *.html
	$(RM) TAGS

realclean: clean
	$(RM) -fr t

TAGS:	$(sources)
	etags $^

$(tdir)/f90.depends: $(sources) $(MODDEP)
	@if [ ! -d $(tdir) ]; then mkdir -p $(tdir); fi
	perl $(MODDEP) --odir $(odir) $(sources) > $(tdir)/f90.depends 
-include $(tdir)/f90.depends
