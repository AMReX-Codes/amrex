#ifndef MAIN_H
#define MAIN_H

#include <AMReX_Gpu.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Particles.H>
#include <AMReX_VisMF.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Vector.H>
#include <cmath>

namespace amrex {

   IntVect cctype(AMREX_D_DECL(0,0,0));
   IntVect ndtype(AMREX_D_DECL(1,1,1));
   IntVect xetype(AMREX_D_DECL(0,1,1));
   IntVect yetype(AMREX_D_DECL(1,0,1));
   IntVect zetype(AMREX_D_DECL(1,1,0));
   IntVect xftype(AMREX_D_DECL(1,0,0));
   IntVect yftype(AMREX_D_DECL(0,1,0));
   IntVect zftype(AMREX_D_DECL(0,0,1));

}

#endif
