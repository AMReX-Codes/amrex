#ifndef Kernels_H_
#define Kernels_H_

// GPU kernels to calculate slopes
#include "slope_K.H"

// GPU kernels to calculate fluxes
#if (AMREX_SPACEDIM == 2)
#include "flux_2d_K.H"
#else
#include "flux_3d_K.H"
#endif

// GPU kernels to update state during advection
#include "Adv_K.H"

// GPU kernels to tag high error cells
#include "tagging_K.H"

// GPU kernel to compute Godunov velocities on each face.
// This is a problem-dependent kernel defined under the Exec directory.
#if (AMREX_SPACEDIM == 2)
#include "face_velocity_2d_K.H"
#else
#include "face_velocity_3d_K.H"
#endif

#endif /* Kernels_H_ */
