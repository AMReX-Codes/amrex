#ifndef _Prob_H_
#define _Prob_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

AMREX_GPU_DEVICE
void initdata(Box const& bx, FArrayBox& phifab,
              GeometryData const & geom);

#endif
