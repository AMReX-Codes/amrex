#
# $Id: GNUmakefile,v 1.1 1998-05-06 17:50:30 lijewski Exp $
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

AmrLib.tex: $(_docHeaders)
	-if `doc++ -h > /dev/null` ; then                                 \
	    doc++ -f -t -j -o AmrLib.tex -ep docxx_squish $(_docHeaders); \
	    sed '/[\]usepackage[{]docxx[}]/d' < AmrLib.tex > .junk;       \
	    mv .junk AmrLib.tex;                                          \
	fi

AmrLib.dvi: AmrLib.tex
	latex AmrLib.tex; latex AmrLib.tex

AmrLib.ps: AmrLib.dvi
	dvips AmrLib.dvi -o AmrLib.ps
