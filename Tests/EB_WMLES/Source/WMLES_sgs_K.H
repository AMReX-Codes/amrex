#ifndef WMLES_SGS_K_H_
#define WMLES_SGS_K_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <AMReX_EBCellFlag.H>
#include <cmath>
#include "WMLES_parm.H"

// ---------------------------------------------------------------------------
// Smagorinsky eddy viscosity: nu_t = (Cs * Delta)^2 * |S|
//   Delta = filter_width_ratio * (dx*dy*dz)^(1/3)
//   |S|   = sqrt(2 * S_ij * S_ij)
// ---------------------------------------------------------------------------

// EB-aware cell-centered first derivative: 2nd-order central difference,
// falling back to a one-sided difference adjacent to covered (cut) cells.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
wmles_grad_eb (int i, int j, int k, int comp, int dir,
               amrex::Array4<amrex::Real const> const& q,
               amrex::Array4<amrex::EBCellFlag const> const& flag,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    using amrex::Real;
    Real qm, q0, qp;
    bool cm, cp;
    if (dir == 0) {
        qm = q(i-1,j,k,comp); q0 = q(i,j,k,comp); qp = q(i+1,j,k,comp);
        cm = flag(i-1,j,k).isCovered(); cp = flag(i+1,j,k).isCovered();
    } else if (dir == 1) {
        qm = q(i,j-1,k,comp); q0 = q(i,j,k,comp); qp = q(i,j+1,k,comp);
        cm = flag(i,j-1,k).isCovered(); cp = flag(i,j+1,k).isCovered();
    } else {
#if (AMREX_SPACEDIM == 3)
        qm = q(i,j,k-1,comp); q0 = q(i,j,k,comp); qp = q(i,j,k+1,comp);
        cm = flag(i,j,k-1).isCovered(); cp = flag(i,j,k+1).isCovered();
#else
        return Real(0.0);
#endif
    }
    if (cm && cp) { return Real(0.0); }
    if (cm)       { return (qp - q0) * dxinv[dir]; }
    if (cp)       { return (q0 - qm) * dxinv[dir]; }
    return (qp - qm) * Real(0.5) * dxinv[dir];
}

// nu_t from the strain-rate magnitude squared (shared by EB / non-EB paths)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
wmles_nu_t_from_S2 (amrex::Real S2,
                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxinv,
                    Parm const& parm) noexcept
{
    using amrex::Real;
    Real strate = std::sqrt(S2);
    Real Delta = parm.sgs_filter_width_ratio;
#if (AMREX_SPACEDIM == 2)
    Delta *= std::sqrt(Real(1.0) / (dxinv[0]*dxinv[1]));
#else
    Delta *= std::cbrt(Real(1.0) / (dxinv[0]*dxinv[1]*dxinv[2]));
#endif
    return (parm.Cs_smag * Delta) * (parm.Cs_smag * Delta) * strate;
}

// EB-safe Smagorinsky nu_t: uses one-sided gradients so covered neighbors
// do not leak garbage values into the strain-rate tensor near cut cells.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
wmles_compute_nu_t (int i, int j, int k,
                    amrex::Array4<amrex::Real const> const& q,
                    amrex::Array4<amrex::EBCellFlag const> const& flag,
                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxinv,
                    Parm const& parm) noexcept
{
    using amrex::Real;

    if (flag(i,j,k).isCovered()) { return Real(0.0); }

    Real dudx = wmles_grad_eb(i,j,k,QU,0,q,flag,dxinv);
    Real dudy = wmles_grad_eb(i,j,k,QU,1,q,flag,dxinv);
    Real dvdx = wmles_grad_eb(i,j,k,QV,0,q,flag,dxinv);
    Real dvdy = wmles_grad_eb(i,j,k,QV,1,q,flag,dxinv);
#if (AMREX_SPACEDIM == 3)
    Real dudz = wmles_grad_eb(i,j,k,QU,2,q,flag,dxinv);
    Real dvdz = wmles_grad_eb(i,j,k,QV,2,q,flag,dxinv);
    Real dwdx = wmles_grad_eb(i,j,k,QW,0,q,flag,dxinv);
    Real dwdy = wmles_grad_eb(i,j,k,QW,1,q,flag,dxinv);
    Real dwdz = wmles_grad_eb(i,j,k,QW,2,q,flag,dxinv);
#endif

    Real S2;
#if (AMREX_SPACEDIM == 2)
    S2 = Real(2.0)*dudx*dudx + Real(2.0)*dvdy*dvdy
       + (dudy+dvdx)*(dudy+dvdx);
#else
    S2 = Real(2.0)*(dudx*dudx + dvdy*dvdy + dwdz*dwdz)
       + (dudy+dvdx)*(dudy+dvdx)
       + (dudz+dwdx)*(dudz+dwdx)
       + (dvdz+dwdy)*(dvdz+dwdy);
#endif

    return wmles_nu_t_from_S2(S2, dxinv, parm);
}

// Non-EB Smagorinsky nu_t: plain 2nd-order central differences.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
wmles_compute_nu_t (int i, int j, int k,
                    amrex::Array4<amrex::Real const> const& q,
                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxinv,
                    Parm const& parm) noexcept
{
    using amrex::Real;

    Real dudx = (q(i+1,j,k,QU) - q(i-1,j,k,QU)) * Real(0.5) * dxinv[0];
    Real dudy = (q(i,j+1,k,QU) - q(i,j-1,k,QU)) * Real(0.5) * dxinv[1];
    Real dvdx = (q(i+1,j,k,QV) - q(i-1,j,k,QV)) * Real(0.5) * dxinv[0];
    Real dvdy = (q(i,j+1,k,QV) - q(i,j-1,k,QV)) * Real(0.5) * dxinv[1];
#if (AMREX_SPACEDIM == 3)
    Real dudz = (q(i,j,k+1,QU) - q(i,j,k-1,QU)) * Real(0.5) * dxinv[2];
    Real dvdz = (q(i,j,k+1,QV) - q(i,j,k-1,QV)) * Real(0.5) * dxinv[2];
    Real dwdx = (q(i+1,j,k,QW) - q(i-1,j,k,QW)) * Real(0.5) * dxinv[0];
    Real dwdy = (q(i,j+1,k,QW) - q(i,j-1,k,QW)) * Real(0.5) * dxinv[1];
    Real dwdz = (q(i,j,k+1,QW) - q(i,j,k-1,QW)) * Real(0.5) * dxinv[2];
#endif

    Real S2;
#if (AMREX_SPACEDIM == 2)
    S2 = Real(2.0)*dudx*dudx + Real(2.0)*dvdy*dvdy
       + (dudy+dvdx)*(dudy+dvdx);
#else
    S2 = Real(2.0)*(dudx*dudx + dvdy*dvdy + dwdz*dwdz)
       + (dudy+dvdx)*(dudy+dvdx)
       + (dudz+dwdx)*(dudz+dwdx)
       + (dvdz+dwdy)*(dvdz+dwdy);
#endif

    return wmles_nu_t_from_S2(S2, dxinv, parm);
}

#endif // WMLES_SGS_K_H_
