
#ifndef WRITEMULTIFAB_H
#define WRITEMULTIFAB_H

//
// $Id: WriteMultiFab.H,v 1.1 1998-03-24 07:06:42 almgren Exp $
//

#include <REAL.H>
#include <Box.H>
#include <BoxList.H>
#include <MultiFab.H>

// WriteMultiFab is a hacked routine to write a MultiFab into a plotfile
// of the same form written by the iamr code, and is to be read
// in with the amrvis tools.  This routine constructs a background
// box to fill the domain, and sets the values to 0.  The MultiFab is
// written as a level=1 data set over the top of the level=0 background

void
WriteMultiFab(char *filename, const MultiFab &phi,
              Real H, const BoxList& bs, const Box &container,
              int ratio=2, Real bgValue=0.0);

void
WriteMultiFab(char *filename, const MultiFab &phi,
              const Real* H, const BoxList& bs, const Box &container,
              int ratio=2, Real bgValue=0.0);

#endif /*WRITEMULTIFAB_H*/
