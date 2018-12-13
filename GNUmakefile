AMREX_HOME  ?= ../amrex
PICSAR_HOME ?= ../picsar
OPENBC_HOME ?= ../openbc_poisson

DEBUG	= FALSE
#DEBUG	= TRUE

#DIM     = 2
DIM	= 3

COMP = gcc
#COMP = intel
#COMP = pgi

TINY_PROFILE   = TRUE
#PROFILE       = TRUE
#COMM_PROFILE  = TRUE
#TRACE_PROFILE = TRUE

STORE_OLD_PARTICLE_ATTRIBS = FALSE

USE_OMP   = TRUE
USE_GPU   = FALSE

EBASE     = main

USE_PYTHON_MAIN = FALSE

USE_SENSEI_INSITU = FALSE

USE_CONDUIT ?= FALSE
USE_ASCENT ?= FALSE

# The following must be absolute paths
CONDUIT_DIR ?= /project/projectdirs/alpine/ghweber/ascent/uberenv_libs/spack/opt/spack/cray-cnl9-haswell/gcc-5.3.0/conduit-master-w6ilqv5kfgt2wgbftzpssjj3lttgkqkv
ASCENT_DIR ?= /project/projectdirs/alpine/ghweber/ascent/uberenv_libs/ascent-install

WarpxBinDir = Bin

USE_PSATD = FALSE

DO_ELECTROSTATIC = FALSE

WARPX_HOME := .
include $(WARPX_HOME)/Source/Make.WarpX
