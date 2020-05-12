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
cns_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state,
              amrex::GeometryData const& geomdata, Parm const& parm, ProbParm const& prob_parm)
{
    using namespace amrex::literals;

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    constexpr amrex::Real pi = 3.14159265358979323846264338327950288_rt;
    const amrex::Real splitx = 0.5_rt*(prob_lo[0]+prob_hi[0]);
    const amrex::Real splity = 0.5_rt*(prob_lo[1]+prob_hi[1]);
    const amrex::Real splitz = 0.5_rt*(prob_lo[2]+prob_hi[2]);
    const amrex::Real L_x = prob_hi[0] - prob_lo[0];
    const amrex::Real presmid = prob_parm.p0_base - prob_parm.rho_1*splitz;

    const amrex::Real z = prob_lo[2] + (k+0.5_rt)*dx[2];
    const amrex::Real y = prob_lo[1] + (j+0.5_rt)*dx[1];
    const amrex::Real x = prob_lo[0] + (i+0.5_rt)*dx[0];

    amrex::Real Pt;
    if (z < splitz) {
        Pt = prob_parm.p0_base - prob_parm.rho_1*z;
    } else {
        Pt = presmid - prob_parm.rho_2*(z-splitz);
    }
    const amrex::Real rhoet = Pt/(parm.eos_gamma-1.0_rt);

    const Real r2d = amrex::min(std::hypot((x-splitx),(y-splity)), 0.5_rt*L_x);
    const Real pertheight = 0.5_rt - 0.01_rt*std::cos(2.0_rt*pi*r2d/L_x);
    const Real rhot = prob_parm.rho_1 + ((prob_parm.rho_2-prob_parm.rho_1)/2.0_rt)*(1.0_rt+std::tanh((z-pertheight)/0.005_rt));

    state(i,j,k,URHO ) = rhot;
    state(i,j,k,UMX  ) = 0.0_rt;
    state(i,j,k,UMY  ) = 0.0_rt;
    state(i,j,k,UMZ  ) = 0.0_rt;
    state(i,j,k,UEINT) = rhoet;
    state(i,j,k,UEDEN) = rhoet;
    state(i,j,k,UTEMP) = 0.0_rt;
}

#endif
