#include <cmath>
#include <limits>

#include <WarpX.H>
#include <WarpXConst.H>
#include <WarpX_f.H>
#ifdef WARPX_USE_PY
#include <WarpX_py.H>
#endif

#ifdef BL_USE_SENSEI_INSITU
#include <AMReX_AmrMeshInSituBridge.H>
#endif

#include <PML_current.H>

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
            WRPX_DAMP_PML(tex.loVect(), tex.hiVect(),
			    tey.loVect(), tey.hiVect(),
			    tez.loVect(), tez.hiVect(),
			    tbx.loVect(), tbx.hiVect(),
    			tby.loVect(), tby.hiVect(),
	    		tbz.loVect(), tbz.hiVect(),
		    	BL_TO_FORTRAN_3D((*pml_E[0])[mfi]),
			    BL_TO_FORTRAN_3D((*pml_E[1])[mfi]),
			    BL_TO_FORTRAN_3D((*pml_E[2])[mfi]),
			    BL_TO_FORTRAN_3D((*pml_B[0])[mfi]),
			    BL_TO_FORTRAN_3D((*pml_B[1])[mfi]),
			    BL_TO_FORTRAN_3D((*pml_B[2])[mfi]),
			    WRPX_PML_TO_FORTRAN(sigba[mfi]));

            if (pml_F) {
                const Box& tnd  = mfi.nodaltilebox();
                WRPX_DAMP_PML_F(tnd.loVect(), tnd.hiVect(),
			        BL_TO_FORTRAN_3D((*pml_F)[mfi]),
			        WRPX_PML_TO_FORTRAN(sigba[mfi]));
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

    BL_PROFILE("WarpX::DampJPML()");

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

            auto const& AMREX_RESTRICT x_lo = sigba[mfi].sigma_cumsum_fac[0].lo();
#if (AMREX_SPACEDIM == 3)
            auto const& AMREX_RESTRICT y_lo = sigba[mfi].sigma_cumsum_fac[1].lo();
            auto const& AMREX_RESTRICT z_lo = sigba[mfi].sigma_cumsum_fac[2].lo();
#else
            int y_lo = 0;
            auto const& AMREX_RESTRICT z_lo = sigba[mfi].sigma_cumsum_fac[1].lo();
#endif

            auto const& AMREX_RESTRICT xs_lo = sigba[mfi].sigma_star_cumsum_fac[0].lo();
#if (AMREX_SPACEDIM == 3)
            auto const& AMREX_RESTRICT ys_lo = sigba[mfi].sigma_star_cumsum_fac[1].lo();
            auto const& AMREX_RESTRICT zs_lo = sigba[mfi].sigma_star_cumsum_fac[2].lo();
#else
            int ys_lo = 0;
            auto const& AMREX_RESTRICT zs_lo = sigba[mfi].sigma_star_cumsum_fac[1].lo();
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

/* \brief Copy the current J from the regular grid to the PML */
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
