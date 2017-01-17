#include <MultiFab.H>
#include <Geometry.H>

#ifndef WARPX_UTIL_H_
#define WARPX_UTIL_H_

void shiftMF(MultiFab& mf, const Geometry& geom, int num_shift, 
	     int dir, const IntVect& nodalflag);

#endif
