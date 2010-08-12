    FC = pathf95
    F90 = pathf95
    CC  = pathcc

    FFLAGS   += -module $(mdir) -I$(mdir) 
    F90FLAGS += -module $(mdir) -I$(mdir)

    ifdef OMP
      F90FLAGS += -mp
      FFLAGS += -mp
      CFLAGS += -mp
    endif

    ifdef USE_HPCTOOLKIT
      HPCLINK = hpclink
      HPCLINK_FLAGS_PATHSCALE = -g
    endif

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
      ifeq ($(findstring nid, $(HOST)), nid)
      #
      # franklin.nersc.gov
      #
        F90FLAGS += -O -fno-second-underscore $(HPCLINK_FLAGS_PATHSCALE)
        FFLAGS   += -O -fno-second-underscore $(HPCLINK_FLAGS_PATHSCALE)
        CFLAGS   += -O -fno-second-underscore $(HPCLINK_FLAGS_PATHSCALE)
      else
        F90FLAGS += -Ofast -fno-second-underscore $(HPCLINK_FLAGS_PATHSCALE)
        FFLAGS   += -Ofast -fno-second-underscore $(HPCLINK_FLAGS_PATHSCALE)
        CFLAGS   += -Ofast -fno-second-underscore $(HPCLINK_FLAGS_PATHSCALE)
      endif
    endif
#   LDFLAGS += -static
    CPPFLAGS += -DBL_HAS_SECOND
