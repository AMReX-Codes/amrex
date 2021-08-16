#ifndef Kernels_H_
#define Kernels_H_

// GPU kernels to tag high error cells
#include "Tagging_nd_K.H"

// GPU kernel to compute Godunov velocities on each face.
// This is a problem-dependent kernel defined under the Exec directory.
#if (AMREX_SPACEDIM == 2)
#include "face_velocity_2d_K.H"
#else
#include "face_velocity_3d_K.H"
#endif

#endif
