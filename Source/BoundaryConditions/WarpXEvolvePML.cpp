/* Copyright 2019 Aurore Blelly, Axel Huebl, Maxence Thevenet
 * Remi Lehe, Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "Utils/WarpXConst.H"
#include "WarpX_PML_kernels.H"
#ifdef WARPX_USE_PY
#   include "Python/WarpX_py.H"
#endif

#include "PML_current.H"

#ifdef BL_USE_SENSEI_INSITU
#   include <AMReX_AmrMeshInSituBridge.H>
#endif

#include <cmath>
#include <limits>


using namespace amrex;

void
WarpX::DampPML ()
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        DampPML(lev);
    }
}

void
WarpX::DampPML (int lev)
{
    DampPML(lev, PatchType::fine);
    if (lev > 0) DampPML(lev, PatchType::coarse);
}

void
WarpX::DampPML (int lev, PatchType patch_type)
{
    if (!do_pml) return;

    WARPX_PROFILE("WarpX::DampPML()");

    if (pml[lev]->ok())
    {
        const auto& pml_E = (patch_type == PatchType::fine) ? pml[lev]->GetE_fp() : pml[lev]->GetE_cp();
        const auto& pml_B = (patch_type == PatchType::fine) ? pml[lev]->GetB_fp() : pml[lev]->GetB_cp();
        const auto& pml_F = (patch_type == PatchType::fine) ? pml[lev]->GetF_fp() : pml[lev]->GetF_cp();
        const auto& sigba = (patch_type == PatchType::fine) ? pml[lev]->GetMultiSigmaBox_fp()
                                                              : pml[lev]->GetMultiSigmaBox_cp();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*pml_E[0], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            const Box& tex  = mfi.tilebox(Ex_nodal_flag);
            const Box& tey  = mfi.tilebox(Ey_nodal_flag);
            const Box& tez  = mfi.tilebox(Ez_nodal_flag);
            const Box& tbx  = mfi.tilebox(Bx_nodal_flag);
            const Box& tby  = mfi.tilebox(By_nodal_flag);
            const Box& tbz  = mfi.tilebox(Bz_nodal_flag);

            auto const& pml_Exfab = pml_E[0]->array(mfi);
            auto const& pml_Eyfab = pml_E[1]->array(mfi);
            auto const& pml_Ezfab = pml_E[2]->array(mfi);
            auto const& pml_Bxfab = pml_B[0]->array(mfi);
            auto const& pml_Byfab = pml_B[1]->array(mfi);
            auto const& pml_Bzfab = pml_B[2]->array(mfi);
            amrex::Real const * AMREX_RESTRICT sigma_fac_x = sigba[mfi].sigma_fac[0].data();
#if (AMREX_SPACEDIM == 3)
            amrex::Real const * AMREX_RESTRICT sigma_fac_y = sigba[mfi].sigma_fac[1].data();
            amrex::Real const * AMREX_RESTRICT sigma_fac_z = sigba[mfi].sigma_fac[2].data();
#else
            amrex::Real const * AMREX_RESTRICT sigma_fac_y = nullptr;
            amrex::Real const * AMREX_RESTRICT sigma_fac_z = sigba[mfi].sigma_fac[1].data();
#endif
            amrex::Real const * AMREX_RESTRICT sigma_star_fac_x = sigba[mfi].sigma_star_fac[0].data();
#if (AMREX_SPACEDIM == 3)
            amrex::Real const * AMREX_RESTRICT sigma_star_fac_y = sigba[mfi].sigma_star_fac[1].data();
            amrex::Real const * AMREX_RESTRICT sigma_star_fac_z = sigba[mfi].sigma_star_fac[2].data();
#else
            amrex::Real const * AMREX_RESTRICT sigma_star_fac_y = nullptr;
            amrex::Real const * AMREX_RESTRICT sigma_star_fac_z = sigba[mfi].sigma_star_fac[1].data();
#endif
            int const x_lo = sigba[mfi].sigma_fac[0].lo();
#if (AMREX_SPACEDIM == 3)
            int const y_lo = sigba[mfi].sigma_fac[1].lo();
            int const z_lo = sigba[mfi].sigma_fac[2].lo();
#else
            int const y_lo = 0;
            int const z_lo = sigba[mfi].sigma_fac[1].lo();
#endif

            amrex::ParallelFor(tex, tey, tez,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                warpx_damp_pml_ex(i,j,k,pml_Exfab,sigma_fac_y,sigma_fac_z,
                                  sigma_star_fac_x,x_lo,y_lo,z_lo);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                warpx_damp_pml_ey(i,j,k,pml_Eyfab,sigma_fac_z,sigma_fac_x,
                                  sigma_star_fac_y,x_lo,y_lo,z_lo);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                warpx_damp_pml_ez(i,j,k,pml_Ezfab,sigma_fac_x,sigma_fac_y,
                                  sigma_star_fac_z,x_lo,y_lo,z_lo);
            });

            amrex::ParallelFor(tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                warpx_damp_pml_bx(i,j,k,pml_Bxfab,sigma_star_fac_y,
                                  sigma_star_fac_z,y_lo,z_lo);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                warpx_damp_pml_by(i,j,k,pml_Byfab,sigma_star_fac_z,
                                  sigma_star_fac_x,z_lo,x_lo);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                warpx_damp_pml_bz(i,j,k,pml_Bzfab,sigma_star_fac_x,
                                  sigma_star_fac_y,x_lo,y_lo);
            });

            if (pml_F) {
               // Note that for warpx_damp_pml_F(), mfi.nodaltilebox is used in
               // the ParallelFor loop and here we use mfi.tilebox.
               /// But, it does not matter because in damp_pml, where
               // nodaltilebox is used, only a simple multiplication is performed.
                const Box& tnd  = mfi.nodaltilebox();
                auto const& pml_F_fab = pml_F->array(mfi);
                amrex::ParallelFor(tnd,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    warpx_damp_pml_F(i,j,k,pml_F_fab,sigma_fac_x,
                                     sigma_fac_y,sigma_fac_z,
                                     x_lo,y_lo,z_lo);
                });

            }
        }
    }
}

void
WarpX::DampJPML ()
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        DampJPML(lev);
    }
}

void
WarpX::DampJPML (int lev)
{
    DampJPML(lev, PatchType::fine);
    if (lev > 0) DampJPML(lev, PatchType::coarse);
}

void
WarpX::DampJPML (int lev, PatchType patch_type)
{
    if (!do_pml) return;
    if (!do_pml_j_damping) return;

    WARPX_PROFILE("WarpX::DampJPML()");

    if (pml[lev]->ok())
    {

        const auto& pml_j = (patch_type == PatchType::fine) ? pml[lev]->Getj_fp() : pml[lev]->Getj_cp();
        const auto& sigba = (patch_type == PatchType::fine) ? pml[lev]->GetMultiSigmaBox_fp()
                                                            : pml[lev]->GetMultiSigmaBox_cp();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*pml_j[0], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            auto const& pml_jxfab = pml_j[0]->array(mfi);
            auto const& pml_jyfab = pml_j[1]->array(mfi);
            auto const& pml_jzfab = pml_j[2]->array(mfi);
            const Real* sigma_cumsum_fac_j_x = sigba[mfi].sigma_cumsum_fac[0].data();
            const Real* sigma_star_cumsum_fac_j_x = sigba[mfi].sigma_star_cumsum_fac[0].data();
#if (AMREX_SPACEDIM == 3)
            const Real* sigma_cumsum_fac_j_y = sigba[mfi].sigma_cumsum_fac[1].data();
            const Real* sigma_star_cumsum_fac_j_y = sigba[mfi].sigma_star_cumsum_fac[1].data();
            const Real* sigma_cumsum_fac_j_z = sigba[mfi].sigma_cumsum_fac[2].data();
            const Real* sigma_star_cumsum_fac_j_z = sigba[mfi].sigma_star_cumsum_fac[2].data();
#else
            const Real* sigma_cumsum_fac_j_y = nullptr;
            const Real* sigma_star_cumsum_fac_j_y = nullptr;
            const Real* sigma_cumsum_fac_j_z = sigba[mfi].sigma_cumsum_fac[1].data();
            const Real* sigma_star_cumsum_fac_j_z = sigba[mfi].sigma_star_cumsum_fac[1].data();
#endif
            const Box& tjx  = mfi.tilebox(jx_nodal_flag);
            const Box& tjy  = mfi.tilebox(jy_nodal_flag);
            const Box& tjz  = mfi.tilebox(jz_nodal_flag);

            int const x_lo = sigba[mfi].sigma_cumsum_fac[0].lo();
#if (AMREX_SPACEDIM == 3)
            int const y_lo = sigba[mfi].sigma_cumsum_fac[1].lo();
            int const z_lo = sigba[mfi].sigma_cumsum_fac[2].lo();
#else
            int const y_lo = 0;
            int const z_lo = sigba[mfi].sigma_cumsum_fac[1].lo();
#endif

            int const xs_lo = sigba[mfi].sigma_star_cumsum_fac[0].lo();
#if (AMREX_SPACEDIM == 3)
            int const ys_lo = sigba[mfi].sigma_star_cumsum_fac[1].lo();
            int const zs_lo = sigba[mfi].sigma_star_cumsum_fac[2].lo();
#else
            int const ys_lo = 0;
            int const zs_lo = sigba[mfi].sigma_star_cumsum_fac[1].lo();
#endif

            amrex::ParallelFor( tjx, tjy, tjz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    damp_jx_pml(i, j, k, pml_jxfab, sigma_star_cumsum_fac_j_x,
                                sigma_cumsum_fac_j_y, sigma_cumsum_fac_j_z,
                                xs_lo,y_lo, z_lo);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    damp_jy_pml(i, j, k, pml_jyfab, sigma_cumsum_fac_j_x,
                                sigma_star_cumsum_fac_j_y, sigma_cumsum_fac_j_z,
                                x_lo,ys_lo, z_lo);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    damp_jz_pml(i, j, k, pml_jzfab, sigma_cumsum_fac_j_x,
                                sigma_cumsum_fac_j_y, sigma_star_cumsum_fac_j_z,
                                x_lo,y_lo, zs_lo);
                }
            );
        }

    }
}

/**
 * \brief Copy the current J from the regular grid to the PML
 */
void
WarpX::CopyJPML ()
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        if (pml[lev]->ok()){
            pml[lev]->CopyJtoPMLs({ current_fp[lev][0].get(),
                                  current_fp[lev][1].get(),
                                  current_fp[lev][2].get() },
                                { current_cp[lev][0].get(),
                                  current_cp[lev][1].get(),
                                  current_cp[lev][2].get() });
        }
    }
}
