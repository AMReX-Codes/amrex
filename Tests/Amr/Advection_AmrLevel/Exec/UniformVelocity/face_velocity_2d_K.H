#ifndef FACE_VELOCITY_2D_H_
#define FACE_VELOCITY_2D_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Array.H>

#include "AmrLevelAdv.H"
#include "Prob_Parm.H"

AMREX_GPU_HOST
AMREX_FORCE_INLINE
void get_face_velocity(const amrex::Real /*time*/,
                       amrex::FArrayBox& vx,
                       amrex::FArrayBox& vy,
                       amrex::GpuArray<amrex::Real,2> /*dx*/,
                       amrex::GpuArray<amrex::Real,2> /*prob_lo*/)
{
    using namespace amrex;

    vx.setVal<RunOn::Device>(AmrLevelAdv::d_prob_parm->adv_vel[0]);
    vy.setVal<RunOn::Device>(AmrLevelAdv::d_prob_parm->adv_vel[1]);

    return;
}

#endif
