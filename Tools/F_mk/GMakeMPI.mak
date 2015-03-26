# MPI reconfiguration...

ifndef MPI
$(error THIS FILE SHOULD ONLY BE LOADED WITH MPI defined)
endif

# Architecture specific changes...
ifeq ($(ARCH),AIX)
  F90 = mpxlf95$(rsuf)
  FC  = mpxlf$(rsuf)
  mpi_libraries = -lmpi
endif

ifeq ($(ARCH),Darwin)

  # specifics for maudib
  ifeq ($(findstring maudib, $(UNAMEN)), maudib)
    FC = mpif90
    F90 = mpif90
    CC = mpicc
    CXX = mpicxx

    MPI_HOME=$(shell dirname `mpicc --showme:libdirs | cut -d" " -f2`)
    LIBRARIES += lmpi_mpifh    

  # attempt at general support for Mac; only works if MPIHOME env var is set
  # and applies some general defaults
  else ifdef MPIHOME
    FC = $(MPIHOME)/bin/mpif90
    F90 = $(MPIHOME)/bin/mpif90
    CC = $(MPIHOME)/bin/mpicc
    CXX = $(MPIHOME)/bin/mpicxx

  # otherwise break
  else
    $(error SORRY, no MPI specification for your Darwin/Mac machine; check BoxLib/Tools/F_mk/GMakeMPI.mak and/or set MPIHOME environment variable)
  endif

  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib

endif

ifeq ($(ARCH),Linux)
  ifeq ($(COMP),g95)
    F90FLAGS += -fno-second-underscore
    FFLAGS   += -fno-second-underscore
    override F_C_LINK := UNDERSCORE
  endif
  ifeq ($(COMP),PathScale)
    FC = mpif90
    F90 = mpif90
    ifdef MPIHOME
      mpi_include_dir = $(MPIHOME)/include
    endif
  endif  # end of PathScale
endif
#
# Try and catch the Lawrencium cluster ...
#
ifeq ($(findstring .scs, $(HOSTNAMEF)), .scs)
    FC  = mpif90
    F90 = mpif90
endif
#
# Host changes ....
#
ifeq ($(findstring intrepid, $(HOSTNAMEF)), intrepid)
    #
    # intrepid.alcf.anl.gov -- we only seem to be able to get the name
    #                          intrepid from hostname -f.  $HOST or
    #                          uname -n don't indicate intrepid
    #

    ifdef OMP
      CC  := mpixlc_r
      FC  := mpixlf95_r -qfixed=72
      F90 := mpixlf95_r
    else  
      CC  := mpixlc
      FC  := mpixlf95 -qfixed=72
      F90 := mpixlf95
    endif

    FFLAGS   += -qmoddir=$(mdir) -I$(mdir)
    F90FLAGS += -qmoddir=$(mdir) -I$(mdir)
    CFLAGS   += -I$(mdir) -Wp,-DBL_AIX

    ifdef NDEBUG
      FFLAGS   += -O2 -qarch=450d -qtune=450
      F90FLAGS += -O2 -qarch=450d -qtune=450
      CFLAGS   += -O2 -qarch=450d -qtune=450
    else
      FFLAGS   += -g -C
      F90FLAGS += -g -C
      CFLAGS   += -g -qcheck=bounds
    endif

    F_C_LINK := LOWERCASE

    # if using the bg* compilers instead of the mpi* wrappers above, you may
    # need these
    #
    #    MPIHOME=/bgsys/drivers/ppcfloor
    #    mpi_include_dir = $(MPIHOME)/arch/include -I$(MPIHOME)/comm/include
    #    mpi_lib_dir = $(MPIHOME)/comm/lib -L$(MPIHOME)/runtime/SPI
    #    mpi_libraries += -lmpich.cnk -ldcmfcoll.cnk -ldcmf.cnk 
    #    mpi_libraries += -lpthread -lrt -lSPI.cna
endif

ifeq ($(findstring surveyor, $(HOSTNAME)), surveyor)
    #
    # surveyor.alcf.anl.gov
    #
    CC  := mpixlc
    FC  := mpixlf95 -qfixed=72
    F90 := mpixlf95

    FFLAGS   := -qmoddir=$(mdir) -I$(mdir)
    F90FLAGS := -qmoddir=$(mdir) -I$(mdir)
    CFLAGS   := -I$(mdir) -Wp,-DBL_AIX

    ifdef NDEBUG
      FFLAGS   += -O2 -qarch=450d -qtune=450
      F90FLAGS += -O2 -qarch=450d -qtune=450
      CFLAGS   += -O2 -qarch=450d -qtune=450
    else
      FFLAGS   += -g -C
      F90FLAGS += -g -C
      CFLAGS   += -g -qcheck=bounds
    endif

    F_C_LINK := LOWERCASE

    # if using the bg* compilers instead of the mpi* wrappers above, you may
    # need these
    #
    #    MPIHOME=/bgsys/drivers/ppcfloor
    #    mpi_include_dir = $(MPIHOME)/arch/include -I$(MPIHOME)/comm/include
    #    mpi_lib_dir = $(MPIHOME)/comm/lib -L$(MPIHOME)/runtime/SPI
    #    mpi_libraries += -lmpich.cnk -ldcmfcoll.cnk -ldcmf.cnk 
    #    mpi_libraries += -lpthread -lrt -lSPI.cna
endif

ifeq ($(findstring cvrsvc, $(HOST)), cvrsvc)
    #
    # carver.nersc.gov
    #
    ifdef MPI
        CXX := mpiCC
        FC  := mpif90
        F90 := mpif90
    endif
endif
ifeq ($(findstring grace, $(HOST)), grace)
    #
    # grace.nersc.gov
    #
    ifdef MPI
        CXX := CC -target=linux
        CC  := cc -target=linux
        FC  := ftn -target=linux
        F90 := ftn -target=linux
    endif
endif
ifeq ($(findstring hopper, $(HOST)), hopper)
    #
    # hopper.nersc.gov
    #
    ifdef MPI
        CXX := CC -target=linux
        CC  := cc -target=linux
        FC  := ftn -target=linux
        F90 := ftn -target=linux
    endif
endif
ifeq ($(findstring edison, $(HOST)), edison)
    #
    # edison.nersc.gov
    #
    ifdef MPI
        CXX := CC
        CC  := cc
        FC  := ftn
        F90 := ftn
    endif
endif
ifeq ($(findstring titan, $(HOST)), titan)
    #
    # titan (Oak Ridge, OLCF machine)
		#
		# Cray machines require you use their compiler wrappers
		# even if you aren't using Cray compiler
    #
    ifdef MPI
        CXX := CC 
        CC  := cc 
        FC  := ftn
        F90 := ftn
    endif
endif
ifeq ($(findstring h2o, $(UNAMEN)), h2o)
    #
    # Blue Waters
    #
    ifdef MPI
        CXX := CC
        FC  := ftn
        F90 := ftn
    endif

    ifeq ($(COMP),Cray)
      FFLAGS += -hpgas_runtime
      F90FLAGS += -hpgas_runtime
      CFLAGS += -hpgas_runtime
    endif
endif


ifeq ($(HOST),cfe3)
  mpi_libraries += -lmpi
endif

ifeq ($(HOST), orga)
  MPIHOME=/usr/local
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_libraries += -lmpich -lmpl -lpthread
  ifeq ($(COMP),g95)
    $(error SORRY NO MPI WITH G95)
  endif
endif

ifeq ($(HOST),naphta)
  MPIHOME=/usr/lib/mpich
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = -L$(MPIHOME)/lib
  mpi_libraries += -lmpich -lpthread
  ifeq ($(COMP),g95)
    $(error SORRY NO MPI WITH G95)
  endif
endif

ifeq ($(HOST),battra.lbl.gov)
  MPIHOME=/usr/lib/mpich
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_libraries += -lmpich -lpthread
endif

ifeq ($(HOST),gigan)
  MPIHOME=/usr/lib/mpich
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_libraries += -lmpich -lpthread
endif

ifeq ($(HOST),gamera)
  MPIHOME=/usr/lib/mpich
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_libraries += -lmpich -lpthread
endif

ifeq ($(HOST),kiryu)
  MPIHOME=/usr/local 
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_libraries += -lmpich -lpthread
  ifeq ($(COMP),g95)
    $(error SORRY NO MPI WITH G95)
  endif

endif

ifeq ($(HOST),manda)
  MPIHOME=/usr/lib/mpich
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_libraries += -lmpich -lpthread
endif
ifeq ($(HOST),hedorah)
  MPIHOME=/usr/lib/mpich
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_libraries += -lmpich -lpthread
endif
ifeq ($(HOST),gojira)
  MPIHOME=/usr/local
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_libraries += -lmpich -lmpichf90 -lpthread
endif
ifeq ($(HOST),atragon)
  MPIHOME=/usr/local
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_libraries += -lmpich -lpthread
endif
ifeq ($(HOST),ebirah)
  MPIHOME=/usr/lib/mpich
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_libraries += -lmpich -lpthread
endif
ifeq ($(HOST),rodan)
  MPIHOME=/usr/lib/mpich
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_libraries += -lmpich -lpthread
endif
ifeq ($(HOST),baragon)
  MPIHOME=/usr/lib/mpich
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_libraries += -lmpich -lpthread
endif
ifeq ($(HOST),megalon)
  MPIHOME=/usr/lib/mpich
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_libraries += -lmpich -lpthread
endif
ifeq ($(HOST),posse)
  MPIHOME=/usr/lib/mpich
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_libraries += -lmpich -lpthread
endif
ifeq ($(HOST),mothra)
  MPIHOME=/usr/lib/mpich
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_libraries += -lmpich -lpthread
endif
ifeq ($(HOST),gimantis)
  MPIHOME=/usr/lib/mpich
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_libraries += -lmpich -lpthread
endif
ifeq ($(HOST),angilas)
  MPIHOME=/usr/local
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_libraries += -lmpich -lmpichf90 -lpthread
endif

ifeq ($(findstring donev, $(HOSTNAME)), donev)
   ifeq ($(MPIVENDOR),OpenMPIv1)
      MPIHOME=/usr/lib64/compat-openmpi
      mpi_include_dir = /usr/include/compat-openmpi-x86_64
      mpi_libraries += -lmpi -lmpi_f77 #-lmpi_f90
      mpi_lib_dir = $(MPIHOME)/lib
   else ifeq ($(MPIVENDOR),OpenMPI) # Latest version
      MPIHOME=$(HOME)/HPC/Libraries/OMPI
      mpi_libraries += -lmpi -lmpi_mpifh # -lmpi -lmpi_f77 #-lmpi_f90
      mpi_include_dir = $(MPIHOME)/include
      mpi_lib_dir = $(MPIHOME)/lib
   else
      MPIHOME=$(HOME)/HPC/Libraries/MPI
      mpi_libraries += -lmpich -lmpichf90 -lpthread
      mpi_include_dir = $(MPIHOME)/include
      mpi_lib_dir = $(MPIHOME)/lib
   endif
else ifeq ($(findstring cims.nyu.edu, $(HOSTNAME)), cims.nyu.edu)
   # OpenMPI v2
   MPIHOME=/usr/lib64/openmpi
   mpi_include_dir = /usr/include/openmpi-x86_64
   mpi_libraries += -lmpi -lmpi_mpifh
   # OpenMPI v1
   #MPIHOME=/usr/lib64/compat-openmpi
   #mpi_include_dir = /usr/include/compat-openmpi-x86_64
   #mpi_libraries += -lmpi -lmpi_f77 #-lmpi_f90
   # Generic stuff:
   mpi_lib_dir = $(MPIHOME)/lib
endif

ifeq ($(HOST),greenstreet)
  MPIHOME=/usr/local/anag/pkg/mpich-1.2.6-intel90
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_include_dir = $(MPIHOME)/include
  mpi_libraries += -lmpich
endif

# generic linux install with MPICH wrappers -- set the 
# BOXLIB_USE_MPI_WRAPPERS environment variable for this
ifdef BOXLIB_USE_MPI_WRAPPERS
    F90 = mpif90
    CXX = mpicxx
endif

ifeq ($(HOST),lookfar)
  MPIHOME=/usr/local
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = -L$(MPIHOME)/lib
  mpi_libraries += -lmpich
  ifeq ($(COMP),Intel)
    FFLAGS += -assume 2underscores
    F90FLAGS += -assume 2underscores
    CFLAGS += -DBL_FORT_USE_DBL_UNDERSCORE
    CFLAGS += -UBL_FORT_USE_UNDERSCORE
  endif
endif

ifeq ($(HOSTNAME),hyades.ucsc.edu)
  F90 := mpiifort
  FC := mpiifort
  fC := mpiifort
endif


ifeq ($(findstring lanl, $(UNAMEN)), lanl)
  F90 := mpif90
  FC := mpif90
  fc := mpif90
  CXX := mpic++
  CC := mpicc
endif
