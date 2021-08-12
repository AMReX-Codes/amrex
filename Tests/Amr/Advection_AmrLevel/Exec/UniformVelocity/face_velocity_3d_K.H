#ifndef FACE_VELOCITY_3D_H_
#define FACE_VELOCITY_3D_H_

#include "Prob_Parm.H"

void get_face_velocity(const amrex::Real /*time*/,
                       amrex::FArrayBox& vx,
                       amrex::FArrayBox& vy,
                       amrex::FArrayBox& vz,
                       const amrex::Real* /*dx*/, const amrex::Real* /*prob_lo*/)
{
    using namespace amrex;

    vx.setVal<RunOn::Device>(AmrLevelAdv::d_prob_parm->adv_vel[0]);
    vy.setVal<RunOn::Device>(AmrLevelAdv::d_prob_parm->adv_vel[1]);
    vz.setVal<RunOn::Device>(AmrLevelAdv::d_prob_parm->adv_vel[2]);

    return;
}

#endif
