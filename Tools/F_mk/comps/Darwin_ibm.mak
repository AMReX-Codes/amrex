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
