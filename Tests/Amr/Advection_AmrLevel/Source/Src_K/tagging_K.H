#ifndef TAGGING_K_H_
#define TAGGING_K_H_

#include <AMReX_Array4.H>
#include <AMReX_Math.H>

// Tag regions of high phi.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
state_error (int i, int j, int k,
             amrex::Array4<char> const& tag,
             amrex::Array4<amrex::Real const> const& state,
             amrex::Real const phi_threshold,
             char const tagval) noexcept
{
    if (state(i,j,k) >= phi_threshold) {
        tag(i,j,k) = tagval;
    }
}

// Tag regions of high phi gradient.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
grad_error (int i, int j, int k,
            amrex::Array4<char> const& tag,
            amrex::Array4<amrex::Real const> const& state,
            amrex::Real const grad_threshold,
            char const tagval) noexcept
{
    amrex::Real ax = amrex::Math::abs(state(i-1,j,k) - state(i,j,k));
    amrex::Real ay = amrex::Math::abs(state(i,j-1,k) - state(i,j,k));
    ax = amrex::max(ax, amrex::Math::abs(state(i,j,k) - state(i+1,j,k)));
    ay = amrex::max(ay, amrex::Math::abs(state(i,j,k) - state(i,j+1,k)));
    amrex::Real az;
    if (AMREX_SPACEDIM == 2) {
        az = 0.0;
    } else {
        az = amrex::Math::abs(state(i,j,k-1) - state(i,j,k));
        az = amrex::max(az, amrex::Math::abs(state(i,j,k) - state(i,j,k+1)));
    }
    if (amrex::max(ax,ay,az) >= grad_threshold) {
        tag(i,j,k) = tagval;
    }
}

#endif /* TAGGING_K_H_ */
