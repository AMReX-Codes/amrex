
#include <cmath>
#include <limits>

#include <WarpX.H>
#include <WarpXConst.H>
#include <WarpX_f.H>
#include <WarpX_K.H>
#include <WarpX_PML_kernels.H>
#ifdef WARPX_USE_PY
#include <WarpX_py.H>
#endif

#ifdef BL_USE_SENSEI_INSITU
#include <AMReX_AmrMeshInSituBridge.H>
#endif

using namespace amrex;

#ifdef WARPX_USE_PSATD

void
WarpX::PushPSATD (amrex::Real a_dt)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dt[lev] == a_dt, "dt must be consistent");
        if (fft_hybrid_mpi_decomposition){
#ifndef AMREX_USE_CUDA // Only available on CPU
            PushPSATD_hybridFFT(lev, a_dt);
#endif
        } else {
            PushPSATD_localFFT(lev, a_dt);
        }
    }
}

void WarpX::PushPSATD_localFFT (int lev, amrex::Real /* dt */)
{
    auto& solver = *spectral_solver_fp[lev];

    // Perform forward Fourier transform
    solver.ForwardTransform(*Efield_fp[lev][0], SpectralFieldIndex::Ex);
    solver.ForwardTransform(*Efield_fp[lev][1], SpectralFieldIndex::Ey);
    solver.ForwardTransform(*Efield_fp[lev][2], SpectralFieldIndex::Ez);
    solver.ForwardTransform(*Bfield_fp[lev][0], SpectralFieldIndex::Bx);
    solver.ForwardTransform(*Bfield_fp[lev][1], SpectralFieldIndex::By);
    solver.ForwardTransform(*Bfield_fp[lev][2], SpectralFieldIndex::Bz);
    solver.ForwardTransform(*current_fp[lev][0], SpectralFieldIndex::Jx);
    solver.ForwardTransform(*current_fp[lev][1], SpectralFieldIndex::Jy);
    solver.ForwardTransform(*current_fp[lev][2], SpectralFieldIndex::Jz);
    solver.ForwardTransform(*rho_fp[lev], SpectralFieldIndex::rho_old, 0);
    solver.ForwardTransform(*rho_fp[lev], SpectralFieldIndex::rho_new, 1);

    // Advance fields in spectral space
    solver.pushSpectralFields();

    // Perform backward Fourier Transform
    solver.BackwardTransform(*Efield_fp[lev][0], SpectralFieldIndex::Ex);
    solver.BackwardTransform(*Efield_fp[lev][1], SpectralFieldIndex::Ey);
    solver.BackwardTransform(*Efield_fp[lev][2], SpectralFieldIndex::Ez);
    solver.BackwardTransform(*Bfield_fp[lev][0], SpectralFieldIndex::Bx);
    solver.BackwardTransform(*Bfield_fp[lev][1], SpectralFieldIndex::By);
    solver.BackwardTransform(*Bfield_fp[lev][2], SpectralFieldIndex::Bz);
}

#endif

void
WarpX::EvolveB (Real a_dt)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        EvolveB(lev, a_dt);
    }
}

void
WarpX::EvolveB (int lev, Real a_dt)
{
    BL_PROFILE("WarpX::EvolveB()");
    EvolveB(lev, PatchType::fine, a_dt);
    if (lev > 0)
    {
        EvolveB(lev, PatchType::coarse, a_dt);
    }
}

void
WarpX::EvolveB (int lev, PatchType patch_type, amrex::Real a_dt)
{
    const int patch_level = (patch_type == PatchType::fine) ? lev : lev-1;
    const std::array<Real,3>& dx = WarpX::CellSize(patch_level);
    Real dtsdx = a_dt/dx[0], dtsdy = a_dt/dx[1], dtsdz = a_dt/dx[2];

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

    // xmin is only used by the picsar kernel with cylindrical geometry,
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

        if (do_nodal) {
            auto const& Bxfab = Bx->array(mfi);
            auto const& Byfab = By->array(mfi);
            auto const& Bzfab = Bz->array(mfi);
            auto const& Exfab = Ex->array(mfi);
            auto const& Eyfab = Ey->array(mfi);
            auto const& Ezfab = Ez->array(mfi);
            amrex::ParallelFor(tbx,
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_bx_nodal(j,k,l,Bxfab,Eyfab,Ezfab,dtsdy,dtsdz);
            });
            amrex::ParallelFor(tby,
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_by_nodal(j,k,l,Byfab,Exfab,Ezfab,dtsdx,dtsdz);
            });
            amrex::ParallelFor(tbz,
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_bz_nodal(j,k,l,Bzfab,Exfab,Eyfab,dtsdx,dtsdy);
            });
        } else {
            // Call picsar routine for each tile
            warpx_push_bvec(
		      tbx.loVect(), tbx.hiVect(),
		      tby.loVect(), tby.hiVect(),
		      tbz.loVect(), tbz.hiVect(),
		      BL_TO_FORTRAN_3D((*Ex)[mfi]),
		      BL_TO_FORTRAN_3D((*Ey)[mfi]),
		      BL_TO_FORTRAN_3D((*Ez)[mfi]),
		      BL_TO_FORTRAN_3D((*Bx)[mfi]),
		      BL_TO_FORTRAN_3D((*By)[mfi]),
		      BL_TO_FORTRAN_3D((*Bz)[mfi]),
                      &dtsdx, &dtsdy, &dtsdz,
                      &xmin, &dx[0],
		      &WarpX::maxwell_fdtd_solver_id);
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

            // adding to move f90 to cpp PML_PORT //
            auto const& pml_Bxfab = pml_B[0]->array(mfi);
            auto const& pml_Byfab = pml_B[1]->array(mfi);
            auto const& pml_Bzfab = pml_B[2]->array(mfi);
            auto const& pml_Exfab = pml_E[0]->array(mfi);
            auto const& pml_Eyfab = pml_E[1]->array(mfi);
            auto const& pml_Ezfab = pml_E[2]->array(mfi);
            if (WarpX::maxwell_fdtd_solver_id == 0) {
               amrex::ParallelFor(tbx,
               [=] AMREX_GPU_DEVICE (int i, int j, int k)
               {
                   warpx_push_pml_bx_yee(i,j,k,pml_Bxfab,pml_Eyfab,pml_Ezfab,
                                        dtsdy,dtsdz);
               });
               amrex::ParallelFor(tby,
               [=] AMREX_GPU_DEVICE (int i, int j, int k)
               {
                   warpx_push_pml_by_yee(i,j,k,pml_Byfab,pml_Exfab,pml_Ezfab,
                                         dtsdx,dtsdz);
               });
               amrex::ParallelFor(tbz,
               [=] AMREX_GPU_DEVICE (int i, int j, int k)
               {
                   warpx_push_pml_bz_yee(i,j,k,pml_Bzfab,pml_Exfab,pml_Eyfab,
                                        dtsdx,dtsdy);
               });
            } else if (WarpX::maxwell_fdtd_solver_id == 1) {
              amrex::Print() << " CKC solver not implemented for CPP run\n";
            }

        }
    }
}

void
WarpX::EvolveE (Real a_dt)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        EvolveE(lev, a_dt);
    }
}

void
WarpX::EvolveE (int lev, Real a_dt)
{
    BL_PROFILE("WarpX::EvolveE()");
    EvolveE(lev, PatchType::fine, a_dt);
    if (lev > 0)
    {
        EvolveE(lev, PatchType::coarse, a_dt);
    }
}

void
WarpX::EvolveE (int lev, PatchType patch_type, amrex::Real a_dt)
{
    const Real mu_c2_dt = (PhysConst::mu0*PhysConst::c*PhysConst::c) * a_dt;
    const Real c2dt = (PhysConst::c*PhysConst::c) * a_dt;

    int patch_level = (patch_type == PatchType::fine) ? lev : lev-1;
    const std::array<Real,3>& dx = WarpX::CellSize(patch_level);
    Real dtsdx_c2 = c2dt/dx[0], dtsdy_c2 = c2dt/dx[1], dtsdz_c2 = c2dt/dx[2];

    MultiFab *Ex, *Ey, *Ez, *Bx, *By, *Bz, *jx, *jy, *jz, *F;
    if (patch_type == PatchType::fine)
    {
        Ex = Efield_fp[lev][0].get();
        Ey = Efield_fp[lev][1].get();
        Ez = Efield_fp[lev][2].get();
        Bx = Bfield_fp[lev][0].get();
        By = Bfield_fp[lev][1].get();
        Bz = Bfield_fp[lev][2].get();
        jx = current_fp[lev][0].get();
        jy = current_fp[lev][1].get();
        jz = current_fp[lev][2].get();
        F  = F_fp[lev].get();
    }
    else if (patch_type == PatchType::coarse)
    {
        Ex = Efield_cp[lev][0].get();
        Ey = Efield_cp[lev][1].get();
        Ez = Efield_cp[lev][2].get();
        Bx = Bfield_cp[lev][0].get();
        By = Bfield_cp[lev][1].get();
        Bz = Bfield_cp[lev][2].get();
        jx = current_cp[lev][0].get();
        jy = current_cp[lev][1].get();
        jz = current_cp[lev][2].get();
        F  = F_cp[lev].get();
    }

    MultiFab* cost = costs[lev].get();
    const IntVect& rr = (lev > 0) ? refRatio(lev-1) : IntVect::TheUnitVector();

    // xmin is only used by the picsar kernel with cylindrical geometry,
    // in which case it is actually rmin.
    const Real xmin = Geom(0).ProbLo(0);

    // Loop through the grids, and over the tiles within each grid
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Ex, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Real wt = amrex::second();

        const Box& tex  = mfi.tilebox(Ex_nodal_flag);
        const Box& tey  = mfi.tilebox(Ey_nodal_flag);
        const Box& tez  = mfi.tilebox(Ez_nodal_flag);

        if (do_nodal) {
            auto const& Exfab = Ex->array(mfi);
            auto const& Eyfab = Ey->array(mfi);
            auto const& Ezfab = Ez->array(mfi);
            auto const& Bxfab = Bx->array(mfi);
            auto const& Byfab = By->array(mfi);
            auto const& Bzfab = Bz->array(mfi);
            auto const& jxfab = jx->array(mfi);
            auto const& jyfab = jy->array(mfi);
            auto const& jzfab = jz->array(mfi);
            amrex::ParallelFor(tex,
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_ex_nodal(j,k,l,Exfab,Byfab,Bzfab,jxfab,mu_c2_dt,dtsdy_c2,dtsdz_c2);
            });
            amrex::ParallelFor(tey,
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_ey_nodal(j,k,l,Eyfab,Bxfab,Bzfab,jyfab,mu_c2_dt,dtsdx_c2,dtsdz_c2);
            });
            amrex::ParallelFor(tez,
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_ez_nodal(j,k,l,Ezfab,Bxfab,Byfab,jzfab,mu_c2_dt,dtsdx_c2,dtsdy_c2);
            });
        } else {
            // Call picsar routine for each tile
            warpx_push_evec(
		      tex.loVect(), tex.hiVect(),
		      tey.loVect(), tey.hiVect(),
		      tez.loVect(), tez.hiVect(),
		      BL_TO_FORTRAN_3D((*Ex)[mfi]),
		      BL_TO_FORTRAN_3D((*Ey)[mfi]),
		      BL_TO_FORTRAN_3D((*Ez)[mfi]),
		      BL_TO_FORTRAN_3D((*Bx)[mfi]),
		      BL_TO_FORTRAN_3D((*By)[mfi]),
		      BL_TO_FORTRAN_3D((*Bz)[mfi]),
		      BL_TO_FORTRAN_3D((*jx)[mfi]),
		      BL_TO_FORTRAN_3D((*jy)[mfi]),
		      BL_TO_FORTRAN_3D((*jz)[mfi]),
		      &mu_c2_dt,
		      &dtsdx_c2, &dtsdy_c2, &dtsdz_c2,
		      &xmin, &dx[0]);
        }

        if (F)
        {
            warpx_push_evec_f(
			  tex.loVect(), tex.hiVect(),
			  tey.loVect(), tey.hiVect(),
			  tez.loVect(), tez.hiVect(),
			  BL_TO_FORTRAN_3D((*Ex)[mfi]),
			  BL_TO_FORTRAN_3D((*Ey)[mfi]),
			  BL_TO_FORTRAN_3D((*Ez)[mfi]),
			  BL_TO_FORTRAN_3D((*F)[mfi]),
                          &dtsdx_c2, &dtsdy_c2, &dtsdz_c2,
			  &WarpX::maxwell_fdtd_solver_id);
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
    }

    if (do_pml && pml[lev]->ok())
    {
        if (F) pml[lev]->ExchangeF(patch_type, F);

        const auto& pml_B = (patch_type == PatchType::fine) ? pml[lev]->GetB_fp() : pml[lev]->GetB_cp();
        const auto& pml_E = (patch_type == PatchType::fine) ? pml[lev]->GetE_fp() : pml[lev]->GetE_cp();
        const auto& pml_F = (patch_type == PatchType::fine) ? pml[lev]->GetF_fp() : pml[lev]->GetF_cp();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*pml_E[0], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            const Box& tex  = mfi.tilebox(Ex_nodal_flag);
            const Box& tey  = mfi.tilebox(Ey_nodal_flag);
            const Box& tez  = mfi.tilebox(Ez_nodal_flag);

            auto const& pml_Exfab = pml_E[0]->array(mfi);
            auto const& pml_Eyfab = pml_E[1]->array(mfi);
            auto const& pml_Ezfab = pml_E[2]->array(mfi);
            auto const& pml_Bxfab = pml_B[0]->array(mfi);
            auto const& pml_Byfab = pml_B[1]->array(mfi);
            auto const& pml_Bzfab = pml_B[2]->array(mfi);

            if (WarpX::maxwell_fdtd_solver_id == 0) {
               amrex::ParallelFor(tex,
               [=] AMREX_GPU_DEVICE (int i, int j, int k)
               {
                   warpx_push_pml_ex_yee(i,j,k,pml_Exfab,pml_Byfab,pml_Bzfab,
                                         dtsdy_c2,dtsdz_c2);
               });
               amrex::ParallelFor(tey,
               [=] AMREX_GPU_DEVICE (int i, int j, int k)
               {
                   warpx_push_pml_ey_yee(i,j,k,pml_Eyfab,pml_Bxfab,pml_Bzfab,
                                         dtsdx_c2,dtsdz_c2);
               });
               amrex::ParallelFor(tez,
               [=] AMREX_GPU_DEVICE (int i, int j, int k)
               {
                   warpx_push_pml_ez_yee(i,j,k,pml_Ezfab,pml_Bxfab,pml_Byfab,
                                         dtsdx_c2,dtsdy_c2);
               });
            }
            else if (WarpX::maxwell_fdtd_solver_id == 1)
            {
              amrex::Print() << " CKC solver not implemented for CPP run\n";
            }

            if (pml_F)
            {

               auto const& pml_F_fab = pml_F->array(mfi);
               amrex::ParallelFor(tex,
               [=] AMREX_GPU_DEVICE (int i, int j, int k)
               {
                   warpx_push_pml_ex_f_yee(i,j,k,pml_Exfab,pml_F_fab,dtsdx_c2);
               });
               amrex::ParallelFor(tey,
               [=] AMREX_GPU_DEVICE (int i, int j, int k)
               {
                   warpx_push_pml_ey_f_yee(i,j,k,pml_Eyfab,pml_F_fab,dtsdy_c2);
               });
               amrex::ParallelFor(tez,
               [=] AMREX_GPU_DEVICE (int i, int j, int k)
               {
                   warpx_push_pml_ez_f_yee(i,j,k,pml_Ezfab,pml_F_fab,dtsdz_c2);
               });
                
            }
        }
    }
}

void
WarpX::EvolveF (Real a_dt, DtType a_dt_type)
{
    if (!do_dive_cleaning) return;

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        EvolveF(lev, a_dt, a_dt_type);
    }
}

void
WarpX::EvolveF (int lev, Real a_dt, DtType a_dt_type)
{
    if (!do_dive_cleaning) return;

    EvolveF(lev, PatchType::fine, a_dt, a_dt_type);
    if (lev > 0) EvolveF(lev, PatchType::coarse, a_dt, a_dt_type);
}

void
WarpX::EvolveF (int lev, PatchType patch_type, Real a_dt, DtType a_dt_type)
{
    if (!do_dive_cleaning) return;

    BL_PROFILE("WarpX::EvolveF()");

    static constexpr Real mu_c2 = PhysConst::mu0*PhysConst::c*PhysConst::c;

    int patch_level = (patch_type == PatchType::fine) ? lev : lev-1;
    const auto& dx = WarpX::CellSize(patch_level);
    const std::array<Real,3> dtsdx {a_dt/dx[0], a_dt/dx[1], a_dt/dx[2]};

    MultiFab *Ex, *Ey, *Ez, *rho, *F;
    if (patch_type == PatchType::fine)
    {
        Ex = Efield_fp[lev][0].get();
        Ey = Efield_fp[lev][1].get();
        Ez = Efield_fp[lev][2].get();
        rho = rho_fp[lev].get();
        F = F_fp[lev].get();
    }
    else
    {
        Ex = Efield_cp[lev][0].get();
        Ey = Efield_cp[lev][1].get();
        Ez = Efield_cp[lev][2].get();
        rho = rho_cp[lev].get();
        F = F_cp[lev].get();
    }

    const int rhocomp = (a_dt_type == DtType::FirstHalf) ? 0 : 1;

    MultiFab src(rho->boxArray(), rho->DistributionMap(), 1, 0);
    ComputeDivE(src, 0, {Ex,Ey,Ez}, dx);
    MultiFab::Saxpy(src, -mu_c2, *rho, rhocomp, 0, 1, 0);
    MultiFab::Saxpy(*F, a_dt, src, 0, 0, 1, 0);

    if (do_pml && pml[lev]->ok())
    {
        const auto& pml_F = (patch_type == PatchType::fine) ? pml[lev]->GetF_fp() : pml[lev]->GetF_cp();
        const auto& pml_E = (patch_type == PatchType::fine) ? pml[lev]->GetE_fp() : pml[lev]->GetE_cp();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*pml_F, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.tilebox();


            auto const& pml_F_fab = pml_F->array(mfi);
            auto const& pml_Exfab = pml_E[0]->array(mfi);
            auto const& pml_Eyfab = pml_E[1]->array(mfi);
            auto const& pml_Ezfab = pml_E[2]->array(mfi);

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                warpx_push_pml_F(i, j, k, pml_F_fab, pml_Exfab,
                                pml_Eyfab, pml_Ezfab,
                                dtsdx[0], dtsdx[1], dtsdx[2]);
            });

        }
    }
}

