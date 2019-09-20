
#include <cmath>
#include <limits>

#include <WarpX.H>
#include <WarpXConst.H>
#include <WarpX_f.H>
#include <WarpX_K.H>
#include <WarpX_FDTD.H>
#ifdef WARPX_USE_PY
#include <WarpX_py.H>
#endif

#include <PML_current.H>

#ifdef BL_USE_SENSEI_INSITU
#include <AMReX_AmrMeshInSituBridge.H>
#endif

using namespace amrex;


void
WarpX::Hybrid_QED_Push (Real a_dt)
{
	if (not do_nodal) {
		try {
			throw do_nodal;
		}
		catch (bool do_nodal) {
			std::cout << "Error: As fo right now the hyrbrid QED algorithm is only compatiable with the nodel lattice scheme.";
			return 0;
		}
	}
	
    for (int lev = 0; lev <= finest_level; ++lev) {
        Hybrid_QED_Push(lev, a_dt);
    }
}

void
WarpX::Hybrid_QED_Push (int lev, Real a_dt)
{
    BL_PROFILE("WarpX::Hybrid_QED_Push()");
    Hybrid_QED_Push(lev, PatchType::fine, a_dt);
    if (lev > 0)
    {
        Hybrid_QED_Push(lev, PatchType::coarse, a_dt);
    }
}

void
WarpX::Hybrid_QED_Push (int lev, PatchType patch_type, amrex::Real a_dt)
{
    const int patch_level = (patch_type == PatchType::fine) ? lev : lev-1;
    const std::array<Real,3>& dx = WarpX::CellSize(patch_level);
    const Real dtsdx = a_dt/dx[0], dtsdy = a_dt/dx[1], dtsdz = a_dt/dx[2];
    const Real dxinv = 1./dx[0];

    MultiFab *Ex, *Ey, *Ez, *Bx, *By, *Bz;
    if (patch_type == PatchType::fine)
    {
        Ex = Efield_fp[lev][0].get();
        Ey = Efield_fp[lev][1].get();
        Ez = Efield_fp[lev][2].get();
        Bx = Bfield_fp[lev][0].get();
        By = Bfield_fp[lev][1].get();
        Bz = Bfield_fp[lev][2].get();
    }
    else
    {
        Ex = Efield_cp[lev][0].get();
        Ey = Efield_cp[lev][1].get();
        Ez = Efield_cp[lev][2].get();
        Bx = Bfield_cp[lev][0].get();
        By = Bfield_cp[lev][1].get();
        Bz = Bfield_cp[lev][2].get();
    }

    MultiFab* cost = costs[lev].get();
    const IntVect& rr = (lev > 0) ? refRatio(lev-1) : IntVect::TheUnitVector();

    // xmin is only used by the kernel for cylindrical geometry,
    // in which case it is actually rmin.
    const Real xmin = Geom(0).ProbLo(0);

    // Loop through the grids, and over the tiles within each grid
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Bx, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Real wt = amrex::second();

        const Box& tbx  = mfi.tilebox(Bx_nodal_flag);
        const Box& tby  = mfi.tilebox(By_nodal_flag);
        const Box& tbz  = mfi.tilebox(Bz_nodal_flag);
		
		const Box& tex  = mfi.tilebox(Ex_nodal_flag);
        const Box& tey  = mfi.tilebox(Ey_nodal_flag);
        const Box& tez  = mfi.tilebox(Ez_nodal_flag);

        auto const& Bxfab = Bx->array(mfi);
        auto const& Byfab = By->array(mfi);
        auto const& Bzfab = Bz->array(mfi);
        auto const& Exfab = Ex->array(mfi);
        auto const& Eyfab = Ey->array(mfi);
        auto const& Ezfab = Ez->array(mfi);
		
		amrex::ParallelFor(tbx, tby, tbz,
		[=] AMREX_GPU_DEVICE (int j, int k, int l)
        {
			std::string direction = "x";
        	warpx_hybrid_QED_push(j,k,l, Exfab, Eyfab, Ezfab, Bxfab, Byfab, Bzfab, dtsdx, dtsdy, dtsdz, direction);
		},
		[=] AMREX_GPU_DEVICE (int j, int k, int l)
		{
			direction = "y";
			warpx_hybrid_QED_push(j,k,l, Exfab, Eyfab, Ezfab, Bxfab, Byfab, Bzfab, dtsdx, dtsdy, dtsdz, direction);
		},
		[=] AMREX_GPU_DEVICE (int j, int k, int l)
		{
			direction = "z";
			warpx_hybrid_QED_push_ez(j,k,l, Exfab, Eyfab, Ezfab, Bxfab, Byfab, Bzfab, dtsdx, dtsdy, dtsdz, direction);
		});
        }

        if (cost) {
            Box cbx = mfi.tilebox(IntVect{AMREX_D_DECL(0,0,0)});
            if (patch_type == PatchType::coarse) cbx.refine(rr);
            wt = (amrex::second() - wt) / cbx.d_numPts();
            auto costfab = cost->array(mfi);
            amrex::ParallelFor(cbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                costfab(i,j,k) += wt;
            });
        }
	

    if (do_pml && pml[lev]->ok())
    {
        const auto& pml_B = (patch_type == PatchType::fine) ? pml[lev]->GetB_fp() : pml[lev]->GetB_cp();
        const auto& pml_E = (patch_type == PatchType::fine) ? pml[lev]->GetE_fp() : pml[lev]->GetE_cp();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*pml_B[0], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            const Box& tbx  = mfi.tilebox(Bx_nodal_flag);
            const Box& tby  = mfi.tilebox(By_nodal_flag);
            const Box& tbz  = mfi.tilebox(Bz_nodal_flag);

            WRPX_PUSH_PML_BVEC(
			     tbx.loVect(), tbx.hiVect(),
			     tby.loVect(), tby.hiVect(),
			     tbz.loVect(), tbz.hiVect(),
			     BL_TO_FORTRAN_3D((*pml_E[0])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_E[1])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_E[2])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_B[0])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_B[1])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_B[2])[mfi]),
                             &dtsdx, &dtsdy, &dtsdz,
			     &WarpX::maxwell_fdtd_solver_id);
        }
    }
}
