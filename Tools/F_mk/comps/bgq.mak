ifeq ($(COMP),IBM)
  # Even though flush() is in the Fortran standard, the IBM XL Fortran
  # compilers for some reason only understand the symbol flush_() (note the
  # underscore). This flag lets us call flush() without the underscore.
  F90FLAGS += -qextname=flush
  FFLAGS   += -qextname=flush

  F90FLAGS += -g
  FFLAGS   += -g
  CFLAGS   += -g
  CXXFLAGS += -g

  # The IBM XL C++ compiler generates the symbol "__IBMCPlusPlusExceptionV3" in
  # C++ files which do not have exception handling. When compiling C++ and
  # Fortran code together, the C++ linker can resolve the symbol, but the
  # Fortran linker cannot. This is why we need this flag in F_mk but not in
  # C_mk.
  CXXFLAGS += -qnoeh

  ifdef OMP
    FC  := bgxlf95_r
    F90 := bgxlf95_r
    CC  := bgxlc_r
    CXX := bgxlC_r

    F90FLAGS += -qsmp=noauto:omp
    FFLAGS   += -qsmp=noauto:omp
    CFLAGS   += -qsmp=noauto:omp
    CXXFLAGS += -qsmp=noauto:omp
  else
    FC  := bgxlf95
    F90 := bgxlf95
    CC  := bgxlc
    CXX := bgxlC
  endif

  F90FLAGS += -I $(mdir) -qmoddir=$(mdir)
  FFLAGS   += -I $(mdir) -qmoddir=$(mdir)

  F_C_LINK := LOWERCASE
  ifdef NDEBUG
    F90FLAGS += -O2
    FFLAGS += -O2
    CFLAGS += -O2
    CXXFLAGS += -O2
  else
    F90FLAGS += -O0
    FFLAGS += -O0
    CFLAGS += -O0
    CXXFLAGS += -O0
  endif

  # manual linkage options suggested by ALCF doc page for mixing F90 with C++
  xtr_libraries += -L/soft/compilers/ibmcmp-nov2012/xlf/bg/14.1/lib64 -lxlopt -lxl -lxlf90 -lxlfmath
  # One of these libraries has "aio_*" functions which libxlf90.a needs.  The
  # MPI wrappers will pick it up but the serial wrappers won't. I just copied
  # this linkage blob from the MPI wrappers by doing "mpixlcxx_r
  # -show".
  xtr_libraries += -L/bgsys/drivers/V1R2M2/ppc64/comm/lib -L/bgsys/drivers/V1R2M2/ppc64/comm/lib -L/bgsys/drivers/V1R2M2/ppc64/comm/lib64 -L/bgsys/drivers/V1R2M2/ppc64/comm/lib -L/bgsys/drivers/V1R2M2/ppc64/spi/lib -L/bgsys/drivers/V1R2M2/ppc64/comm/sys/lib -L/bgsys/drivers/V1R2M2/ppc64/spi/lib -L/bgsys/drivers/V1R2M2/ppc64/comm/sys/lib -L/bgsys/drivers/V1R2M2/ppc64/comm/lib64 -L/bgsys/drivers/V1R2M2/ppc64/comm/lib -L/bgsys/drivers/V1R2M2/ppc64/spi/lib -L/bgsys/drivers/V1R2M2/ppc64/comm/lib-lopa-xl -lmpl-xl -lpami-gcc -lSPI -lSPI_cnk -lrt -lpthread -lstdc++ -lpthread
endif
