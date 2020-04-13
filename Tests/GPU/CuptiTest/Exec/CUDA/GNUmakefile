# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../../../

DEBUG     = FALSE
USE_MPI   = TRUE 
USE_OMP   = FALSE
USE_CUDA  = TRUE
USE_CUPTI = TRUE
COMP      = gcc
DIM       = 3 

TINY_PROFILE   = TRUE

Bpack   := ./Make.package
Blocs   := .

include ../../Make.CUPTI
