#ifndef _WritePlotFile_H_
#define _WritePlotFile_H_

#include <aString.H>
#include <MultiFab.H>

// Write a MultiFab to a plotfile named "name" using the HyperCLawV1.1 parallel
// plotfile data format
void
writePlotFile (const char*     name,
	       const MultiFab& mf,
	       const Geometry& geom,
	       const IntVect&  refRatio,
	       Real            bgVal = 0.0);

#endif

