  ifdef OMP
    FC  := xlf95_r
    F90 := xlf95_r
    CC  := xlc_r
    CXX := xlC_r

    F90FLAGS += -qsmp=noauto:omp
    FFLAGS   += -qsmp=noauto:omp
    CFLAGS   += -qsmp=noauto:omp
    CXXFLAGS += -qsmp=noauto:omp
  else
    FC  := xlf95
    F90 := xlf95
    CC  := xlc
    CXX := xlC
  endif 

  F90FLAGS += -I $(mdir) -qmoddir=$(mdir)
  FFLAGS   += -I $(mdir) -qmoddir=$(mdir)

  F_C_LINK := LOWERCASE
  ifdef NDEBUG
    F90FLAGS += -O 
    FFLAGS += -O 
    CFLAGS += -O
    CXXFLAGS += -O
  else
    F90FLAGS += -g 
    FFLAGS += -g 
    CFLAGS += -g
    CXXFLAGS += -g
  endif

ifeq ($(UNAMEN),fen)
 HOST = nyblue
 F_C_LINK := LOWERCASE

 # to run something on the front end nodes (fen's) we can't cross compile
 # to specify that we want to compile for the fen, set MACHINE=fen at runtime
 ifeq ($(MACHINE), fen)

   CC  := /opt/ibmcmp/vac/bg/9.0/bin/xlc
   CXX := /opt/ibmcmp/vac/bg/9.0/bin/xlC
   FC  := /opt/ibmcmp/xlf/bg/11.1/bin/xlf
   F90 := /opt/ibmcmp/xlf/bg/11.1/bin/xlf95

   F90FLAGS += -qmoddir=$(mdir) -I$(mdir)
   F90FLAGS += -O2 -qarch=440d -qtune=440

   FFLAGS += -qmoddir=$(mdir) -I$(mdir) -qfixed=72
   FFLAGS += -O2 -qarch=440d -qtune=440

   CFLAGS += -O2 -qarch=440d -qtune=440
   CFLAGS += -I$(mdir) -Wp,-DBL_AIX

   CXXFLAGS += -O2 -qarch=440d -qtune=440
   CXXFLAGS += -I$(mdir) -Wp,-DBL_AIX

 else

   CC  := mpixlc
   CXX := mpixlC
   FC  := mpixlf90
   F90 := mpixlf95

   F90FLAGS += -qmoddir=$(mdir) -I$(mdir)
   F90FLAGS += -O2 -qarch=440d -qtune=440

   FFLAGS += -qmoddir=$(mdir) -I$(mdir) -qfixed=72
   FFLAGS += -O2 -qarch=440d -qtune=440

   CFLAGS += -O2 -qarch=440d -qtune=440
   CFLAGS += -I$(mdir) -Wp,-DBL_AIX

   CXXFLAGS += -O2 -qarch=440d -qtune=440
   CXXFLAGS += -I$(mdir) -Wp,-DBL_AIX
 endif

endif

ifeq ($(UNAMEN),fenp)
 HOST = nyblue
 F_C_LINK := LOWERCASE
 CC  := mpixlc
 CXX := mpixlC
 FC  := mpixlf95
 F90 := mpixlf95

 F90FLAGS += -qmoddir=$(mdir) -I$(mdir)
 F90FLAGS += -O2 -qarch=450d -qtune=450

 FFLAGS += -qmoddir=$(mdir) -I$(mdir) -qfixed=72
 FFLAGS += -O2 -qarch=450d -qtune=450

 CFLAGS += -O2 -qarch=450d -qtune=450
 CFLAGS += -I$(mdir) -Wp,-DBL_AIX

 CXXFLAGS += -O2 -qarch=450d -qtune=450
 CXXFLAGS += -I$(mdir) -Wp,-DBL_AIX

endif
