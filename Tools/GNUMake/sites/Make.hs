#
# For NEC Aurora
#

ifeq ($(USE_MPI),TRUE)
  CC  = ncc
  CXX = nc++
  FC  = nfort
  F90 = nfort
  LIBRARIES += -lmpichf90
endif
