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
