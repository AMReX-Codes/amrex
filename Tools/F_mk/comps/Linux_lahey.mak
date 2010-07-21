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
