#
# $Id: GNUmakefile,v 1.1 1998-05-06 17:52:43 lijewski Exp $
#

#
# Temp stuff for doing html docs.
#
_docHeaders := $(sort $(wildcard *.H))

_htmlDir := html

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

BndryLib.tex: $(_docHeaders)
	-if `doc++ -h > /dev/null` ; then                                   \
	    doc++ -f -t -j -o BndryLib.tex -ep docxx_squish $(_docHeaders); \
	    sed '/[\]usepackage[{]docxx[}]/d' < BndryLib.tex > .junk;       \
	    mv .junk BndryLib.tex;                                          \
	fi

BndryLib.dvi: BndryLib.tex
	latex BndryLib.tex; latex BndryLib.tex

BndryLib.ps: BndryLib.dvi
	dvips BndryLib.dvi -o BndryLib.ps
