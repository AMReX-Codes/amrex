AMREX_HOME ?= ../../../../../../../..

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

#DIM        = 2
DIM       = 3
CPPFLAGS += -DUSE_PERILLA -DMPICH_SKIP_MPICXX -std=c++11  -DUSE_PERILLA_PTHREADS
#-DUSE_PERILLA_OMP
#-DUSE_PERILLA_PTHREADS 
#-DPERILLA_USE_UPCXX

COMP	   = intel
FCOMP	   = gfortran

USE_MPI    = TRUE
USE_OMP    = FALSE

Bpack   := ./Make.package 
Blocs   := . 

include ../Make.Adv.mpi.pthreads
#include ../Make.Adv.mpi.omp
#include ../Make.Adv.upcxx.pthreads
