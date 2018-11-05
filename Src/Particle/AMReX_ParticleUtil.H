#ifndef AMREX_PARTICLEUTIL_H_
#define AMREX_PARTICLEUTIL_H_

#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_Gpu.H>

namespace amrex
{
  AMREX_GPU_HOST_DEVICE
  int getTileIndex (const IntVect& iv, const Box& box, const bool a_do_tiling, 
		    const IntVect& a_tile_size, Box& tbx);  
}

#endif // include guard
