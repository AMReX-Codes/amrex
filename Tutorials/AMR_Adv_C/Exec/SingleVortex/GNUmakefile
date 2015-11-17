BOXLIB_HOME ?= /path/to/BoxLib
SDCLIB_HOME ?= /path/to/SDCLib

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

DIM        = 2

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI       = FALSE
USE_REACT     = TRUE
USE_DIFFUSION = TRUE
USE_SDCLIB    = FALSE

Bpack   := ./Make.package 
Blocs   := . 

ifeq ($(USE_SDCLIB), TRUE)
Bpack += $(SDCLIB_HOME)/Make.package  
Blocs += $(SDCLIB_HOME)/include $(SDCLIB_HOME)/src
endif

include ../Make.ADR
