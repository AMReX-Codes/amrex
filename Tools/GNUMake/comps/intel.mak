CXX = icpx
has_intel_llvm_compiler := $(shell $(CXX) -dumpversion)

ifdef has_intel_llvm_compiler
  $(info Using intel-llvm)
  include $(AMREX_HOME)/Tools/GNUMake/comps/intel-llvm.mak
else
  $(info Using intel-classic)
  include $(AMREX_HOME)/Tools/GNUMake/comps/intel-classic.mak
endif
