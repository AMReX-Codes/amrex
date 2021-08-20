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

    state(i,j,k,UMX  ) = Real(0.0);
    state(i,j,k,UMY  ) = Real(0.0);
    state(i,j,k,UMZ  ) = Real(0.0);

    state(i,j,k,UTEMP) = prob_parm.interior_T;
    state(i,j,k,UEINT) = prob_parm.interior_p / (parm.eos_gamma - Real(1.0));
    state(i,j,k,UEDEN) = state(i,j,k,UEINT);

    state(i,j,k,URHO ) = state(i,j,k,UEINT) / (parm.cv * prob_parm.interior_T);
}

#endif
