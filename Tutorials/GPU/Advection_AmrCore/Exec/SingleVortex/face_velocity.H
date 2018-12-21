#ifndef face_velocity_H_
#define face_velocity_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>

using namespace amrex;

    AMREX_GPU_DEVICE
    void get_face_velocity_psi(Box const& bx,
                               const Real time,
                               FArrayBox& psifab,
                               GeometryData const& geom);

    AMREX_GPU_DEVICE
    void get_face_velocity_x(Box const& bx,
                             FArrayBox& xvelfab,
                             const FArrayBox& psifab,
                             GeometryData const& geom);

    AMREX_GPU_DEVICE
    void get_face_velocity_y(Box const& bx,
                             FArrayBox& yvelfab,
                             const FArrayBox& psifab,
                             GeometryData const& geom);

    AMREX_GPU_DEVICE
    void get_face_velocity_z(Box const& bx,
                             FArrayBox& zvelfab,
                             const FArrayBox& psifab,
                             GeometryData const& geom);

#endif
