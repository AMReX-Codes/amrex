#ifndef FACE_VELOCITY_2D_H_
#define FACE_VELOCITY_2D_H_

#include "Prob_Parm.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity(const amrex::Real /*time*/,
                       amrex::FArrayBox& vx,
                       amrex::FArrayBox& vy,
                       const amrex::Real* /*dx*/, const amrex::Real* /*prob_lo*/)
{
    using namespace amrex;

    vx.setVal<RunOn::Device>(AmrLevelAdv::d_prob_parm->adv_vel[0]);
    vy.setVal<RunOn::Device>(AmrLevelAdv::d_prob_parm->adv_vel[1]);

    return;
}

#endif
