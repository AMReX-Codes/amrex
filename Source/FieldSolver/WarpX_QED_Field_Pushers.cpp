/* Copyright 2019-2020 Glenn Richardson, Maxence Thevenet, Revathi Jambunathan, Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "Utils/WarpXConst.H"
#include "WarpX_QED_K.H"
#include "BoundaryConditions/WarpX_PML_kernels.H"
#include "BoundaryConditions/PML_current.H"
#include "WarpX_FDTD.H"

#ifdef WARPX_USE_PY
#   include "Python/WarpX_py.H"
#endif

#ifdef BL_USE_SENSEI_INSITU
#   include <AMReX_AmrMeshInSituBridge.H>
#endif

#include <cmath>
#include <limits>


using namespace amrex;


void
WarpX::Hybrid_QED_Push (amrex::Vector<amrex::Real> a_dt)
{
    if (WarpX::do_nodal == 0) {
        Print()<<"The do_nodal flag is tripped.\n";
        try{
            throw "Error: The Hybrid QED method is currently only compatible with the nodal scheme.\n";
        }
        catch (const char* msg) {
            std::cerr << msg << std::endl;
            exit(0);
        }
    }
    for (int lev = 0; lev <= finest_level; ++lev) {
        Hybrid_QED_Push(lev, a_dt[lev]);
    }
}

void
WarpX::Hybrid_QED_Push (int lev, Real a_dt)
{
    WARPX_PROFILE("WarpX::Hybrid_QED_Push()");
    Hybrid_QED_Push(lev, PatchType::fine, a_dt);
    if (lev > 0)
    {
        Hybrid_QED_Push(lev, PatchType::coarse, a_dt);
    }
}

void
WarpX::Hybrid_QED_Push (int lev, PatchType patch_type, Real a_dt)
{
    const int patch_level = (patch_type == PatchType::fine) ? lev : lev-1;
    const std::array<Real,3>& dx_vec= WarpX::CellSize(patch_level);
    const Real dx = dx_vec[0];
    const Real dy = dx_vec[1];
    const Real dz = dx_vec[2];

    MultiFab *Ex, *Ey, *Ez, *Bx, *By, *Bz, *Jx, *Jy, *Jz;
    if (patch_type == PatchType::fine)
    {
        Ex = Efield_fp[lev][0].get();
        Ey = Efield_fp[lev][1].get();
        Ez = Efield_fp[lev][2].get();
        Bx = Bfield_fp[lev][0].get();
        By = Bfield_fp[lev][1].get();
        Bz = Bfield_fp[lev][2].get();
        Jx = current_fp[lev][0].get();
        Jy = current_fp[lev][1].get();
        Jz = current_fp[lev][2].get();
    }
    else
    {
        Ex = Efield_cp[lev][0].get();
        Ey = Efield_cp[lev][1].get();
        Ez = Efield_cp[lev][2].get();
        Bx = Bfield_cp[lev][0].get();
        By = Bfield_cp[lev][1].get();
        Bz = Bfield_cp[lev][2].get();
        Jx = current_cp[lev][0].get();
        Jy = current_cp[lev][1].get();
        Jz = current_cp[lev][2].get();
    }

    amrex::Vector<amrex::Real>* cost = WarpX::getCosts(lev);

    // Loop through the grids, and over the tiles within each grid
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Bx, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
        }
        Real wt = amrex::second();

        // Get boxes for E, B, and J

        const Box& tbx = mfi.tilebox(Bx->ixType().toIntVect());
        const Box& tby = mfi.tilebox(By->ixType().toIntVect());
        const Box& tbz = mfi.tilebox(Bz->ixType().toIntVect());

        const Box& tex = mfi.tilebox(Ex->ixType().toIntVect());
        const Box& tey = mfi.tilebox(Ey->ixType().toIntVect());
        const Box& tez = mfi.tilebox(Ez->ixType().toIntVect());

        const Box& tjx = mfi.tilebox(Jx->ixType().toIntVect());
        const Box& tjy = mfi.tilebox(Jy->ixType().toIntVect());
        const Box& tjz = mfi.tilebox(Jz->ixType().toIntVect());

        // Get field arrays
        auto const& Bxfab = Bx->array(mfi);
        auto const& Byfab = By->array(mfi);
        auto const& Bzfab = Bz->array(mfi);
        auto const& Exfab = Ex->array(mfi);
        auto const& Eyfab = Ey->array(mfi);
        auto const& Ezfab = Ez->array(mfi);
        auto const& Jxfab = Jx->array(mfi);
        auto const& Jyfab = Jy->array(mfi);
        auto const& Jzfab = Jz->array(mfi);

        // Define grown box with 1 ghost cell for finite difference stencil
        const Box& gex = amrex::grow(tex,1);
        const Box& gey = amrex::grow(tey,1);
        const Box& gez = amrex::grow(tez,1);

        // Temporary arrays for electric field, protected by Elixir on GPU
        FArrayBox tmpEx_fab(gex,1);
        Elixir tmpEx_eli = tmpEx_fab.elixir();
        auto const& tmpEx = tmpEx_fab.array();

        FArrayBox tmpEy_fab(gey,1);
        Elixir tmpEy_eli = tmpEy_fab.elixir();
        auto const& tmpEy = tmpEy_fab.array();

        FArrayBox tmpEz_fab(gez,1);
        Elixir tmpEz_eli = tmpEz_fab.elixir();
        auto const& tmpEz = tmpEz_fab.array();

        // Copy electric field to temporary arrays
        AMREX_PARALLEL_FOR_4D(
            gex, 1, i, j, k, n,
            { tmpEx(i,j,k,n) = Exfab(i,j,k,n); }
        );

        AMREX_PARALLEL_FOR_4D(
            gey, 1, i, j, k, n,
            { tmpEy(i,j,k,n) = Eyfab(i,j,k,n); }
        );

        AMREX_PARALLEL_FOR_4D(
            gez, 1, i, j, k, n,
            { tmpEz(i,j,k,n) = Ezfab(i,j,k,n); }
        );

        // Make local copy of xi, to use on device.
        const Real xi_c2 = WarpX::quantum_xi_c2;

        // Apply QED correction to electric field, using temporary arrays.
        amrex::ParallelFor(
            tbx,
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_hybrid_QED_push(j,k,l, Exfab, Eyfab, Ezfab, Bxfab, Byfab,
                    Bzfab, tmpEx, tmpEy, tmpEz, Jxfab, Jyfab, Jzfab, dx, dy, dz,
                    a_dt, xi_c2);
            }
        );

        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
            wt = amrex::second() - wt;
            amrex::HostDevice::Atomic::Add( &(*cost)[mfi.index()], wt);
        }
    }
}
