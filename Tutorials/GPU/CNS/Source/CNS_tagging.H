#ifndef CNS_TAGGING_H_
#define CNS_TAGGING_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_TagBox.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE
inline
void
cns_tag_denerror (int i, int j, int k,
                  amrex::Array4<char> const& tag,
                  amrex::Array4<amrex::Real const> const& rho,
                  amrex::Real dengrad_threshold, char tagval)
{
    amrex::Real ax = std::abs(rho(i+1,j,k) - rho(i,j,k));
    amrex::Real ay = std::abs(rho(i,j+1,k) - rho(i,j,k));
    amrex::Real az = std::abs(rho(i,j,k+1) - rho(i,j,k));
    ax = amrex::max(ax,std::abs(rho(i,j,k) - rho(i-1,j,k)));
    ay = amrex::max(ay,std::abs(rho(i,j,k) - rho(i,j-1,k)));
    az = amrex::max(az,std::abs(rho(i,j,k) - rho(i,j,k-1)));
    if (amrex::max(ax,ay,az) >= dengrad_threshold) {
        tag(i,j,k) = tagval;
    }
}

#endif
