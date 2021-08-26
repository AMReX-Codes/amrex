#ifndef FACE_VELOCITY_3D_H_
#define FACE_VELOCITY_3D_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Array.H>

#include "AmrLevelAdv.H"
#include "Prob_Parm.H"

AMREX_GPU_HOST
AMREX_FORCE_INLINE
void get_face_velocity(const amrex::Real /*time*/,
                       amrex::FArrayBox& vx,
                       amrex::FArrayBox& vy,
                       amrex::FArrayBox& vz,
                       amrex::GpuArray<amrex::Real,3> /*dx*/,
                       amrex::GpuArray<amrex::Real,3> /*prob_lo*/)
{
    using namespace amrex;

    vx.setVal<RunOn::Device>(AmrLevelAdv::d_prob_parm->adv_vel[0]);
    vy.setVal<RunOn::Device>(AmrLevelAdv::d_prob_parm->adv_vel[1]);
    vz.setVal<RunOn::Device>(AmrLevelAdv::d_prob_parm->adv_vel[2]);

    return;
}

#endif
