NDEBUG := t
MPI    :=
OMP    :=

MKVERBOSE :=t

COMP := gfortran

# some routines need an eos/network (i.e. to compute thermodynamic
# quantities.  If that is the case, set NEED_EOS_NETWORK := t
NEED_EOS_NETWORK :=

# define the location of the MAESTRO root directory
MAESTRO_TOP_DIR := $(MAESTRO_HOME)


# include the main Makefile stuff
AMREX_HOME := ../../..
include $(AMREX_HOME)/Tools/F_mk/GMakedefs.mak

# core BoxLib directories
BOXLIB_CORE := Src/F_BaseLib

# other packages needed for data_processing
Fmdirs :=


# directories containing files that are 'include'-d via Fortran
Fmincludes :=

ifdef NEED_EOS_NETWORK
  EOS_TOP_DIR := $(MAESTRO_TOP_DIR)/Microphysics/EOS
  NETWORK_TOP_DIR := $(MAESTRO_TOP_DIR)/Microphysics/networks

  EOS_DIR := helmeos
  NETWORK_DIR := ignition_simple

  MICROPHYS_CORE := $(MAESTRO_TOP_DIR)/Microphysics/EOS

  MICROPHYS_CORE += $(EOS_TOP_DIR)/$(EOS_DIR) \
                    $(NETWORK_TOP_DIR)/$(NETWORK_DIR)

  ifeq ($(findstring helmeos, $(EOS_DIR)), helmeos)
    Fmincludes += Microphysics/helmeos
  endif

  f90sources += probin.f90

  include $(NETWORK_TOP_DIR)/$(strip $(NETWORK_DIR))/NETWORK_REQUIRES

  ifdef NEED_VODE
    Fmdirs += Util/VODE
  endif

  ifdef NEED_BLAS
    Fmdirs += Util/BLAS
  endif

endif


# any MAESTRO stuff
Fmpack := $(foreach dir, $(Fmdirs), $(MAESTRO_TOP_DIR)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(Fmdirs), $(MAESTRO_TOP_DIR)/$(dir))

# BoxLib stuff
Fmpack += $(foreach dir, $(BOXLIB_CORE), $(AMREX_HOME)/$(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(BOXLIB_CORE), $(AMREX_HOME)/$(dir))

# Microphysics
Fmpack += $(foreach dir, $(MICROPHYS_CORE), $(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(MICROPHYS_CORE), $(dir))

# any include files
Fmincs := $(foreach dir, $(Fmincludes), $(MAESTRO_TOP_DIR)/$(dir))



# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files
VPATH_LOCATIONS += $(Fmlocs)

# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)

#programs += faverage
#programs += plt_compare_diff_grids
programs += fcompare
#programs += fextract
#programs += fextrema
#programs += fIDLdump
#programs += fsnapshot2d
#programs += fsnapshot3d
#programs += ftime
#programs += fintgvar2d
#programs += fdumpdata2d
#programs += fcoarsen
#programs += fwritecontents2d

ifeq ($(DO_ALL), TRUE)
  ptmp := $(wildcard f*.f90)
  programs := $(basename $(ptmp))
endif

all: $(pnames)

# probin stuff
ifdef NEED_EOS_NETWORK
  PROBIN_PARAMETER_DIRS :=
  EXTERN_PARAMETER_DIRS += $(MICROPHYS_CORE)
  EXTERN_PARAMETERS := $(shell $(AMREX_HOME)/Tools/F_scripts/findparams.py $(EXTERN_PARAMETER_DIRS))

  PROBIN_TEMPLATE := $(MAESTRO_TOP_DIR)/Util/parameters/dummy.probin.template
endif

probin.f90: $(PROBIN_PARAMETERS) $(EXTERN_PARAMETERS) $(PROBIN_TEMPLATE)
	@echo " "
	$(AMREX_HOME)/Tools/F_scripts/write_probin.py \
           -t $(PROBIN_TEMPLATE) -o probin.f90 -n probin \
           --pa "$(PROBIN_PARAMETERS)" --pb "$(EXTERN_PARAMETERS)"


include $(AMREX_HOME)/Tools/F_mk/GMakerules.mak


# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)
