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
    using amrex::Real;

    const Real* prob_lo = geomdata.ProbLo();
    const Real* prob_hi = geomdata.ProbHi();
    const Real* dx      = geomdata.CellSize();

    constexpr Real pi = Real(3.14159265358979323846264338327950288);
    const Real splitx = Real(0.5)*(prob_lo[0]+prob_hi[0]);
    const Real splity = Real(0.5)*(prob_lo[1]+prob_hi[1]);
    const Real splitz = Real(0.5)*(prob_lo[2]+prob_hi[2]);
    const Real L_x = prob_hi[0] - prob_lo[0];
    const Real presmid = prob_parm.p0_base - prob_parm.rho_1*splitz;

    const Real z = prob_lo[2] + (k+Real(0.5))*dx[2];
    const Real y = prob_lo[1] + (j+Real(0.5))*dx[1];
    const Real x = prob_lo[0] + (i+Real(0.5))*dx[0];

    Real Pt;
    if (z < splitz) {
        Pt = prob_parm.p0_base - prob_parm.rho_1*z;
    } else {
        Pt = presmid - prob_parm.rho_2*(z-splitz);
    }
    const Real rhoet = Pt/(parm.eos_gamma-Real(1.0));

    const Real r2d = amrex::min(std::hypot((x-splitx),(y-splity)), Real(0.5)*L_x);
    const Real pertheight = Real(0.5) - Real(0.01)*std::cos(Real(2.0)*pi*r2d/L_x);
    const Real rhot = prob_parm.rho_1 + ((prob_parm.rho_2-prob_parm.rho_1)/Real(2.0))*(Real(1.0)+std::tanh((z-pertheight)/Real(0.005)));

    state(i,j,k,URHO ) = rhot;
    state(i,j,k,UMX  ) = Real(0.0);
    state(i,j,k,UMY  ) = Real(0.0);
    state(i,j,k,UMZ  ) = Real(0.0);
    state(i,j,k,UEINT) = rhoet;
    state(i,j,k,UEDEN) = rhoet;
    state(i,j,k,UTEMP) = Real(0.0);
}

#endif
