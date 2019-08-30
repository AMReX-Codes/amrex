#include <cmath>
#include <limits>

#include <WarpX.H>
#include <WarpXConst.H>
#include <WarpX_f.H>
#include <WarpX_PML_kernels.H>
#ifdef WARPX_USE_PY
#include <WarpX_py.H>
#endif

#ifdef BL_USE_SENSEI_INSITU
#include <AMReX_AmrMeshInSituBridge.H>
#endif

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

    BL_PROFILE("WarpX::DampPML()");

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
            auto const& AMREX_RESTRICT sigma_fac_x = sigba[mfi].sigma_fac[0].data();
#if (AMREX_SPACEDIM == 3) 
            auto const& AMREX_RESTRICT sigma_fac_y = sigba[mfi].sigma_fac[1].data();
            auto const& AMREX_RESTRICT sigma_fac_z = sigba[mfi].sigma_fac[2].data();
#else
            Real* AMREX_RESTRICT sigma_fac_y = nullptr;
            auto const& AMREX_RESTRICT sigma_fac_z = sigba[mfi].sigma_fac[1].data();
#endif
            auto const& AMREX_RESTRICT sigma_star_fac_x = sigba[mfi].sigma_star_fac[0].data();
#if (AMREX_SPACEDIM == 3) 
            auto const& AMREX_RESTRICT sigma_star_fac_y = sigba[mfi].sigma_star_fac[1].data();
            auto const& AMREX_RESTRICT sigma_star_fac_z = sigba[mfi].sigma_star_fac[2].data();
#else
            Real* AMREX_RESTRICT sigma_star_fac_y = nullptr;
            auto const& AMREX_RESTRICT sigma_star_fac_z = sigba[mfi].sigma_star_fac[1].data();
#endif
            auto const& AMREX_RESTRICT x_lo = sigba[mfi].sigma_fac[0].lo();
#if (AMREX_SPACEDIM == 3) 
            auto const& AMREX_RESTRICT y_lo = sigba[mfi].sigma_fac[1].lo();
            auto const& AMREX_RESTRICT z_lo = sigba[mfi].sigma_fac[2].lo();
#else
            int y_lo = 0;
            auto const& AMREX_RESTRICT z_lo = sigba[mfi].sigma_fac[1].lo();
#endif
 
            amrex::ParallelFor(tex,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                warpx_damp_pml_ex(i,j,k,pml_Exfab,sigma_fac_y,sigma_fac_z,
                                  sigma_star_fac_x,x_lo,y_lo,z_lo);
            });

            amrex::ParallelFor(tey,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                warpx_damp_pml_ey(i,j,k,pml_Eyfab,sigma_fac_z,sigma_fac_x,
                                  sigma_star_fac_y,x_lo,y_lo,z_lo);
            });

            amrex::ParallelFor(tez,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                warpx_damp_pml_ez(i,j,k,pml_Ezfab,sigma_fac_x,sigma_fac_y,
                                  sigma_star_fac_z,x_lo,y_lo,z_lo);
            });

            amrex::ParallelFor(tbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                warpx_damp_pml_bx(i,j,k,pml_Bxfab,sigma_star_fac_y,
                                  sigma_star_fac_z,y_lo,z_lo);
            });

            amrex::ParallelFor(tby,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                warpx_damp_pml_by(i,j,k,pml_Byfab,sigma_star_fac_z,
                                  sigma_star_fac_x,z_lo,x_lo);
            });

            amrex::ParallelFor(tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                warpx_damp_pml_bz(i,j,k,pml_Bzfab,sigma_star_fac_x,
                                  sigma_star_fac_y,x_lo,y_lo);
            });


            if (pml_F) {
               // Note that for warpx_damp_pml_F(), mfi.nodaltilebox is used in the ParallelFor loop and here we use mfi.tilebox. But, it does not matter because in damp_pml, where nodaltilebox is used, only a simple multiplication is performed. //   
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
