#
# $Id: GNUmakefile,v 1.25 2010-10-01 21:54:57 marc Exp $
#

BOXLIB_HOME = ../..

TOP = $(BOXLIB_HOME)

PRECISION = DOUBLE
DEBUG     = TRUE
DIM       = 3
COMP      = 
USE_MPI   = TRUE
#NAMESPACE = TRUE
NAMESPACE = FALSE

#
# The base name of the library we're building.
#
LBASE = box
#
# Where libraries and include files will be installed.
#
INSTALL_ROOT = $(BOXLIB_HOME)

ifeq ($(BL_NO_FORT),TRUE)
  DEFINES += -DBL_NO_FORT
endif

ifeq ($(BL_MANGLE_SYMBOLS_WITH_DIM),TRUE)
  ifeq ($(DIM),3)
    DEFINES += -DBox=Box3D
    DEFINES += -DBoxArray=BoxArray3D
    DEFINES += -DBoxDomain=BoxDomain3D
    DEFINES += -DBoxList=BoxList3D
    DEFINES += -DDistributionMapping=DistributionMapping3D
    DEFINES += -DFABio=FABio3D
    DEFINES += -DFABio_8bit=FABio_8bit3D
    DEFINES += -DFABio_ascii=FABio_ascii3D
    DEFINES += -DFABio_binary=FABio_binary3D
    DEFINES += -DFArrayBox=FArrayBox3D
    DEFINES += -DFabArrayBase=FabArrayBase3D
    DEFINES += -DFabArrayId=FabArrayId3D
    DEFINES += -DFillBoxId=FillBoxId3D
    DEFINES += -DIndexType=IndexType3D
    DEFINES += -DIntVect=IntVect3D
    DEFINES += -DMFIter=MFIter3D
    DEFINES += -DMultiFab=MultiFab3D
    DEFINES += -DMultiFabCopyDescriptor=MultiFabCopyDescriptor3D
    DEFINES += -DMultiFabId=MultiFabId3D
    DEFINES += -DParallelDescriptor=ParallelDescriptor3D
    DEFINES += -DParmParse=ParmParse3D
    DEFINES += -DRealBox=RealBox3D
    DEFINES += -DVisMF=VisMF3D
  else
    DEFINES += -DBox=Box2D
    DEFINES += -DBoxArray=BoxArray2D
    DEFINES += -DBoxDomain=BoxDomain2D
    DEFINES += -DBoxList=BoxList2D
    DEFINES += -DDistributionMapping=DistributionMapping2D
    DEFINES += -DFABio=FABio2D
    DEFINES += -DFABio_8bit=FABio_8bit2D
    DEFINES += -DFABio_ascii=FABio_ascii2D
    DEFINES += -DFABio_binary=FABio_binary2D
    DEFINES += -DFArrayBox=FArrayBox2D
    DEFINES += -DFabArrayBase=FabArrayBase2D
    DEFINES += -DFabArrayId=FabArrayId2D
    DEFINES += -DFillBoxId=FillBoxId2D
    DEFINES += -DIndexType=IndexType2D
    DEFINES += -DIntVect=IntVect2D
    DEFINES += -DMFIter=MFIter2D
    DEFINES += -DMultiFab=MultiFab2D
    DEFINES += -DMultiFabCopyDescriptor=MultiFabCopyDescriptor2D
    DEFINES += -DMultiFabId=MultiFabId2D
    DEFINES += -DParallelDescriptor=ParallelDescriptor2D
    DEFINES += -DParmParse=ParmParse2D
    DEFINES += -DRealBox=RealBox2D
    DEFINES += -DVisMF=VisMF2D
  endif
endif

include $(TOP)/Tools/C_mk/Make.defs ./Make.package

ifeq ($(KCC_VERSION),3.3)
CXXFLAGS += --one_instantiation_per_object
endif

all: $(optionsLib)
#
# `install' does NOT try to build and install the HTML stuff.
# doc++ isn't widely distributed enough to assume that.
# For now you must invoke `install_html' separately.
#
install: all

include $(TOP)/Tools/C_mk/Make.rules

#
# Temp stuff for doing html docs.
#
_docHeaders := $(sort $(wildcard *.H))

_htmlDir := html

_installDir = $(instHtmlDir)/boxlib
#
# This unconditionally installs stuff from ./html into _installDir.
#
install_html: html
	-if [ ! -d $(_installDir) ]; then \
              mkdir -p $(_installDir);    \
              chmod 755 $(_installDir);   \
        fi
	-cd $(_installDir); $(RM) *
	if [ -d $(_htmlDir) ]; then                     \
	    cd $(_htmlDir);                             \
	    for f in *; do                              \
	    	    $(INSTALL_DATA) $$f $(_installDir); \
	    done;                                       \
	fi
#
# This will make the HTML files in the directory html provided
# any of the relevant .H files have changed more recently than the
# directory itself.
#
html: $(_docHeaders)
	-if `doc++ -h > /dev/null` ; then                                 \
		$(RM) -r $(_htmlDir);                                     \
		mkdir -p $(_htmlDir);                                     \
		doc++ -f -j -B Banner.html -d $(_htmlDir) $(_docHeaders); \
	fi

BoxLib.tex: $(_docHeaders)
	-if `doc++ -h > /dev/null` ; then                                 \
	    doc++ -f -t -j -o BoxLib.tex -ep docxx_squish $(_docHeaders); \
	    sed '/[\]usepackage[{]docxx[}]/d' < BoxLib.tex > .junk;       \
	    mv .junk BoxLib.tex;                                          \
	fi

BoxLib.dvi: BoxLib.tex
	latex BoxLib.tex; latex BoxLib.tex

BoxLib.ps: BoxLib.dvi
	dvips BoxLib.dvi -o BoxLib.ps
#
# Install the PostScript file into ./ps
#
install_ps: BoxLib.ps
	-if [ ! -d $(INSTALL_ROOT)/ps ]; then \
              mkdir -p $(INSTALL_ROOT)/ps;    \
              chmod 755 $(INSTALL_ROOT)/ps;   \
        fi
	$(INSTALL_DATA) BoxLib.ps $(INSTALL_ROOT)/ps

WRAPPERGEN=$(BOXLIB_HOME)/util/wrappergen/wrappergen
BLMpi.cpp: BLMpi.w
	$(WRAPPERGEN) -c++ -o BLMpi.cpp -w BLMpi.w

clean::
	$(RM) -r html
