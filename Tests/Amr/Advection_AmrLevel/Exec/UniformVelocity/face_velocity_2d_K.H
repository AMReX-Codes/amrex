#ifndef FACE_VELOCITY_2D_H_
#define FACE_VELOCITY_2D_H_

#include "probdata.H"

void get_face_velocity(const amrex::Real /*time*/,
                       amrex::FArrayBox& vx,
                       amrex::FArrayBox& vy,
                       const amrex::Real* /*dx*/, const amrex::Real* /*prob_lo*/)
{
    vx.setVal(adv_vel[0]);
    vy.setVal(adv_vel[1]);

    return;
}

#endif
