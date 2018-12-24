#ifndef CNS_TAGGING_H_
#define CNS_TAGGING_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_TagBox.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE
inline
void
cns_tag_denerror (amrex::Box const& bx, amrex::TagBox& tagfab, amrex::FArrayBox const& rhofab,
                  amrex::Real dengrad_threshold, char tagval)
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto tag = tagfab.view(lo);
    const auto rho = rhofab.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
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
        }
    }
}

#endif
