#
# Setup for using clang/flang
#
CXX = armclang++
CC  = armclang
FC  = armflang
F90 = armflang

CXXFLAGS =
CFLAGS   =
FFLAGS   =
F90FLAGS =

########################################################################

#armclang_version       = $(shell $(CXX) --version | head -1 | sed -e 's/.*version *\([0-9]\+\.[0-9]\+\).*/\1/')
#armclang_major_version = $(shell $(CXX) --version | head -1 | sed -e 's/.*version *\([0-9]\+\.[0-9]\+\).*/\1/' | sed -e 's/\.[0-9]\+//')
#armclang_minor_version = $(shell $(CXX) --version | head -1 | sed -e 's/.*version.*\([0-9]\+\.[0-9]\+\).*/\1/' | sed -e 's/.*\.//')

########################################################################

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -g -O0 -ftrapv
  CFLAGS   += -g -O0 -ftrapv
  FFLAGS   += -g -O0 -ftrapv
  F90FLAGS += -g -O0 -ftrapv

else

  CXXFLAGS += -g -O3
  CFLAGS   += -g -O3
  FFLAGS   += -g -O3
  F90FLAGS += -g -O3

endif

ifeq ($(WARN_ALL),TRUE)
  warning_flags = -Wall -Wextra -Wno-sign-compare -Wunreachable-code -Wnull-dereference
  warning_flags += -Wfloat-conversion -Wextra-semi

  ifneq ($(USE_CUDA),TRUE)
    warning_flags += -Wpedantic
  endif

  ifneq ($(WARN_SHADOW),FALSE)
    warning_flags += -Wshadow
  endif

  CXXFLAGS += $(warning_flags) -Woverloaded-virtual
  CFLAGS += $(warning_flags)
endif

ifeq ($(WARN_ERROR),TRUE)
  CXXFLAGS += -Werror
  CFLAGS += -Werror
endif

# disable some warnings
CXXFLAGS += -Wno-c++17-extensions

########################################################################

ifdef CXXSTD
  CXXSTD := $(strip $(CXXSTD))
else
  CXXSTD := c++17
endif

CXXFLAGS += -std=$(CXXSTD)
CFLAGS   += -std=c11

FMODULES = -J$(fmoddir) -I $(fmoddir)

########################################################################

GENERIC_COMP_FLAGS =

ifeq ($(EXPORT_DYNAMIC),TRUE)
  CPPFLAGS += -DAMREX_EXPORT_DYNAMIC
  LIBRARIES += -Xlinker -export_dynamic
  GENERIC_COMP_FLAGS += -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer
endif

ifeq ($(THREAD_SANITIZER),TRUE)
  GENERIC_COMP_FLAGS += -fsanitize=thread
endif
ifeq ($(FSANITIZER),TRUE)
  GENERIC_COMP_FLAGS += -fsanitize=address -fsanitize=undefined
endif

ifeq ($(USE_OMP),TRUE)
  GENERIC_COMP_FLAGS += -fopenmp
endif

CXXFLAGS += $(GENERIC_COMP_FLAGS) -pthread
CFLAGS   += $(GENERIC_COMP_FLAGS)
FFLAGS   += $(GENERIC_COMP_FLAGS)
F90FLAGS += $(GENERIC_COMP_FLAGS)

########################################################################

ifneq ($(BL_NO_FORT),TRUE)
  override XTRALIBS += -lflang
endif

ifeq ($(FSANITIZER),TRUE)
  override XTRALIBS += -lubsan
endif
