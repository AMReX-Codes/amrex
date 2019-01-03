#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include "CNS_index_macros.H"
#include "CNS_parm.H"
#include "cns_prob_parm.H"

AMREX_GPU_DEVICE
inline
void
cns_initdata (amrex::Box const& bx, amrex::FArrayBox& statefab, amrex::GeometryData const& geomdata)
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto state = statefab.view(lo);
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    constexpr amrex::Real pi = 3.14159265358979323846264338327950288;
    const amrex::Real splitx = 0.5*(prob_lo[0]+prob_hi[0]);
    const amrex::Real splity = 0.5*(prob_lo[1]+prob_hi[1]);
    const amrex::Real splitz = 0.5*(prob_lo[2]+prob_hi[2]);
    const amrex::Real L_x = prob_hi[0] - prob_lo[0];
    const amrex::Real presmid = ProbParm::p0_base - ProbParm::rho_1*splitz;

    for         (int k = 0; k < len.z; ++k) {
        const amrex::Real z = prob_lo[2] + (k+lo.z+0.5)*dx[2];
        for     (int j = 0; j < len.y; ++j) {
            const amrex::Real y = prob_lo[1] + (j+lo.y+0.5)*dx[1];
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                const amrex::Real x = prob_lo[0] + (i+lo.x+0.5)*dx[0];

                amrex::Real Pt;
                if (z < splitz) {
                    Pt = ProbParm::p0_base - ProbParm::rho_1*z;
                } else {
                    Pt = presmid - ProbParm::rho_2*(z-splitz);
                }
                const amrex::Real rhoet = Pt/(Parm::eos_gamma-1.0);

                const Real r2d = amrex::min(std::hypot((x-splitx),(y-splity)), 0.5*L_x);
                const Real pertheight = 0.5 - 0.01*std::cos(2.0*pi*r2d/L_x);
                const Real rhot = ProbParm::rho_1 + ((ProbParm::rho_2-ProbParm::rho_1)/2.0)*(1.0+std::tanh((z-pertheight)/0.005));

                state(i,j,k,URHO ) = rhot;
                state(i,j,k,UMX  ) = 0.0;
                state(i,j,k,UMY  ) = 0.0;
                state(i,j,k,UMZ  ) = 0.0;
                state(i,j,k,UEINT) = rhoet;
                state(i,j,k,UEDEN) = rhoet;
                state(i,j,k,UTEMP) = 0.0;
            }
        }
    }
}

#endif
