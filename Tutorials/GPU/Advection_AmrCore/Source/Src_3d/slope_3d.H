#ifndef slope_3d_H_
#define slope_3d_H_

#include <AMReX_BLFort.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

using namespace amrex;

AMREX_GPU_DEVICE
void slopex2(Box const& bx,
             const FArrayBox &qfab,
             FArrayBox &dqfab);

AMREX_GPU_DEVICE
void slopex4(Box const& bx,
             const FArrayBox &qfab,
             const FArrayBox &dqfab,
             FArrayBox &dq4fab);

AMREX_GPU_DEVICE
void slopey2(Box const& bx,
             const FArrayBox &qfab,
             FArrayBox &dqfab);

AMREX_GPU_DEVICE
void slopey4(Box const& bx,
             const FArrayBox &qfab,
             const FArrayBox &dqfab,
             FArrayBox &dq4fab);

AMREX_GPU_DEVICE
void slopez2(Box const& bx,
             const FArrayBox &qfab,
             FArrayBox &dqfab);

AMREX_GPU_DEVICE
void slopez4(Box const& bx,
             const FArrayBox &qfab,
             const FArrayBox &dqfab,
             FArrayBox &dq4fab);

#endif
