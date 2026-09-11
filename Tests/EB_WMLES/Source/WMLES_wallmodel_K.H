#ifndef WMLES_WALLMODEL_K_H_
#define WMLES_WALLMODEL_K_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <AMReX_EBCellFlag.H>
#include "WMLES_parm.H"

// ===========================================================================
// Wall-model kernel for WMLES (ODE-based, Phase 4 implementation target)
//
// For each cut-cell adjacent to an EB wall:
//   1. Sample velocity along wall-normal direction at N sample heights
//   2. Solve 1D ODE: d/dn[(nu + nu_t) du/dn] = 0
//   3. Return wall shear stress tau_w
//
// CURRENT STATUS: Placeholder — returns zero stress.
// Phase 4 will implement: sample-and-solve ODE with Thomas tridiagonal.
// ===========================================================================

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
wmles_wall_shear_stress (int i, int j, int k,
                         amrex::Array4<amrex::Real const> const& q,
                         amrex::Array4<amrex::Real const> const& bndrynorm,
                         amrex::Array4<amrex::Real const> const& bndrycent,
                         amrex::Array4<amrex::Real const> const& bndryarea,
                         amrex::Array4<amrex::EBCellFlag const> const& flag,
                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxinv,
                         amrex::Real nu_mol,
                         Parm const& parm) noexcept
{
    using amrex::Real;

    amrex::GpuArray<Real,AMREX_SPACEDIM> tau_w = {};

    // Only compute on cut-cells with EB boundary
    if (!flag(i,j,k).isSingleValued()) { return tau_w; }
    if (bndryarea(i,j,k) <= Real(0.0)) { return tau_w; }

    // TODO Phase 4: Sample velocity along wall-normal direction
    // auto const n = GpuArray<Real,3>{bndrynorm(i,j,k,0), bndrynorm(i,j,k,1), bndrynorm(i,j,k,2)};
    // auto const bc = GpuArray<Real,3>{bndrycent(i,j,k,0), bndrycent(i,j,k,1), bndrycent(i,j,k,2)};
    // Sample u_|| at h_j = boundary_centroid + (j * dy) * normal
    // Solve ODE d/dn[(nu+nu_t) du/dn] = 0 → tau_w

    return tau_w;
}

#endif // WMLES_WALLMODEL_K_H_
