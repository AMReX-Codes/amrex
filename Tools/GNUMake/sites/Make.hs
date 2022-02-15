#
# For Hamburg Observatory computers
#
HS_MACHINES := aurora

ifneq ($(which_computer), $(findstring $(which_computer), $(HS_MACHINES)))
  $(error Unknown HS computer, $(which_computer))
endif

ifeq ($(which_computer),$(filter $(which_computer),aurora))

  ifeq ($(USE_MPI),TRUE)
    CC  = mpincc
    CXX = mpinc++
    FC  = mpinfort
    F90 = mpinfort
  endif

endif