
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

#ifdef WARPX_USE_PSATD
namespace {
    void
    PushPSATDSinglePatch (
        SpectralSolver& solver,
        std::array<std::unique_ptr<amrex::MultiFab>,3>& Efield,
        std::array<std::unique_ptr<amrex::MultiFab>,3>& Bfield,
        std::array<std::unique_ptr<amrex::MultiFab>,3>& current,
        std::unique_ptr<amrex::MultiFab>& rho ) {

        using Idx = SpectralFieldIndex;

        // Perform forward Fourier transform
        solver.ForwardTransform(*Efield[0], Idx::Ex);
        solver.ForwardTransform(*Efield[1], Idx::Ey);
        solver.ForwardTransform(*Efield[2], Idx::Ez);
        solver.ForwardTransform(*Bfield[0], Idx::Bx);
        solver.ForwardTransform(*Bfield[1], Idx::By);
        solver.ForwardTransform(*Bfield[2], Idx::Bz);
        solver.ForwardTransform(*current[0], Idx::Jx);
        solver.ForwardTransform(*current[1], Idx::Jy);
        solver.ForwardTransform(*current[2], Idx::Jz);
        solver.ForwardTransform(*rho, Idx::rho_old, 0);
        solver.ForwardTransform(*rho, Idx::rho_new, 1);
        // Advance fields in spectral space
        solver.pushSpectralFields();
        // Perform backward Fourier Transform
        solver.BackwardTransform(*Efield[0], Idx::Ex);
        solver.BackwardTransform(*Efield[1], Idx::Ey);
        solver.BackwardTransform(*Efield[2], Idx::Ez);
        solver.BackwardTransform(*Bfield[0], Idx::Bx);
        solver.BackwardTransform(*Bfield[1], Idx::By);
        solver.BackwardTransform(*Bfield[2], Idx::Bz);
    }
}

void
WarpX::PushPSATD (amrex::Real a_dt)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dt[lev] == a_dt, "dt must be consistent");
        if (fft_hybrid_mpi_decomposition){
#ifdef WARPX_USE_PSATD_HYBRID
            PushPSATD_hybridFFT(lev, a_dt);
#endif
        } else {
            PushPSATD_localFFT(lev, a_dt);
        }

        // Evolve the fields in the PML boxes
        if (do_pml && pml[lev]->ok()) {
            pml[lev]->PushPSATD();
        }
    }
}

void
WarpX::PushPSATD_localFFT (int lev, amrex::Real /* dt */)
{
    // Update the fields on the fine and coarse patch
    PushPSATDSinglePatch( *spectral_solver_fp[lev],
        Efield_fp[lev], Bfield_fp[lev], current_fp[lev], rho_fp[lev] );
    if (spectral_solver_cp[lev]) {
        PushPSATDSinglePatch( *spectral_solver_cp[lev],
             Efield_cp[lev], Bfield_cp[lev], current_cp[lev], rho_cp[lev] );
    }
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

        auto const& Bxfab = Bx->array(mfi);
        auto const& Byfab = By->array(mfi);
        auto const& Bzfab = Bz->array(mfi);
        auto const& Exfab = Ex->array(mfi);
        auto const& Eyfab = Ey->array(mfi);
        auto const& Ezfab = Ez->array(mfi);
        if (do_nodal) {
            amrex::ParallelFor(tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_bx_nodal(j,k,l,Bxfab,Eyfab,Ezfab,dtsdy,dtsdz);
            },
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_by_nodal(j,k,l,Byfab,Exfab,Ezfab,dtsdx,dtsdz);
            },
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_bz_nodal(j,k,l,Bzfab,Exfab,Eyfab,dtsdx,dtsdy);
            });
        } else if (WarpX::maxwell_fdtd_solver_id == 0) {
            amrex::ParallelFor(tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_bx_yee(j,k,l,Bxfab,Eyfab,Ezfab,dtsdy,dtsdz);
            },
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_by_yee(j,k,l,Byfab,Exfab,Ezfab,dtsdx,dtsdz);
            },
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_bz_yee(j,k,l,Bzfab,Exfab,Eyfab,dtsdx,dtsdy,dxinv,xmin);
            });
        } else if (WarpX::maxwell_fdtd_solver_id == 1) {
            Real betaxy, betaxz, betayx, betayz, betazx, betazy;
            Real gammax, gammay, gammaz;
            Real alphax, alphay, alphaz;
            warpx_calculate_ckc_coefficients(dtsdx, dtsdy, dtsdz,
                                             betaxy, betaxz, betayx, betayz, betazx, betazy,
                                             gammax, gammay, gammaz,
                                             alphax, alphay, alphaz);
            amrex::ParallelFor(tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_bx_ckc(j,k,l,Bxfab,Eyfab,Ezfab,
                                  betaxy, betaxz, betayx, betayz, betazx, betazy,
                                  gammax, gammay, gammaz,
                                  alphax, alphay, alphaz);
            },
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_by_ckc(j,k,l,Byfab,Exfab,Ezfab,
                                  betaxy, betaxz, betayx, betayz, betazx, betazy,
                                  gammax, gammay, gammaz,
                                  alphax, alphay, alphaz);
            },
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_bz_ckc(j,k,l,Bzfab,Exfab,Eyfab,
                                  betaxy, betaxz, betayx, betayz, betazx, betazy,
                                  gammax, gammay, gammaz,
                                  alphax, alphay, alphaz);
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

    const int patch_level = (patch_type == PatchType::fine) ? lev : lev-1;
    const std::array<Real,3>& dx = WarpX::CellSize(patch_level);
    const Real dtsdx_c2 = c2dt/dx[0], dtsdy_c2 = c2dt/dx[1], dtsdz_c2 = c2dt/dx[2];
    const Real dxinv = 1./dx[0];

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

    // xmin is only used by the kernel for cylindrical geometry,
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

        auto const& Exfab = Ex->array(mfi);
        auto const& Eyfab = Ey->array(mfi);
        auto const& Ezfab = Ez->array(mfi);
        auto const& Bxfab = Bx->array(mfi);
        auto const& Byfab = By->array(mfi);
        auto const& Bzfab = Bz->array(mfi);
        auto const& jxfab = jx->array(mfi);
        auto const& jyfab = jy->array(mfi);
        auto const& jzfab = jz->array(mfi);

        if (do_nodal) {
            amrex::ParallelFor(tex, tey, tez,
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_ex_nodal(j,k,l,Exfab,Byfab,Bzfab,jxfab,mu_c2_dt,dtsdy_c2,dtsdz_c2);
            },
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_ey_nodal(j,k,l,Eyfab,Bxfab,Bzfab,jyfab,mu_c2_dt,dtsdx_c2,dtsdz_c2);
            },
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_ez_nodal(j,k,l,Ezfab,Bxfab,Byfab,jzfab,mu_c2_dt,dtsdx_c2,dtsdy_c2);
            });
        } else {
            amrex::ParallelFor(tex, tey, tez,
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_ex_yee(j,k,l,Exfab,Byfab,Bzfab,jxfab,mu_c2_dt,dtsdy_c2,dtsdz_c2);
            },
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_ey_yee(j,k,l,Eyfab,Bxfab,Bzfab,jyfab,mu_c2_dt,dtsdx_c2,dtsdz_c2,xmin);
            },
            [=] AMREX_GPU_DEVICE (int j, int k, int l)
            {
                warpx_push_ez_yee(j,k,l,Ezfab,Bxfab,Byfab,jzfab,mu_c2_dt,dtsdx_c2,dtsdy_c2,dxinv,xmin);
            });
        }

        if (F)
        {
            auto const& Ffab = F->array(mfi);
            if (WarpX::maxwell_fdtd_solver_id == 0) {
                amrex::ParallelFor(tex, tey, tez,
                [=] AMREX_GPU_DEVICE (int j, int k, int l)
                {
                    warpx_push_ex_f_yee(j,k,l,Exfab,Ffab,dtsdx_c2);
                },
                [=] AMREX_GPU_DEVICE (int j, int k, int l)
                {
                    warpx_push_ey_f_yee(j,k,l,Eyfab,Ffab,dtsdy_c2);
                },
                [=] AMREX_GPU_DEVICE (int j, int k, int l)
                {
                    warpx_push_ez_f_yee(j,k,l,Ezfab,Ffab,dtsdz_c2);
                });
            }
            else if (WarpX::maxwell_fdtd_solver_id == 1) {
                Real betaxy, betaxz, betayx, betayz, betazx, betazy;
                Real gammax, gammay, gammaz;
                Real alphax, alphay, alphaz;
                warpx_calculate_ckc_coefficients(dtsdx_c2, dtsdy_c2, dtsdz_c2,
                                                 betaxy, betaxz, betayx, betayz, betazx, betazy,
                                                 gammax, gammay, gammaz,
                                                 alphax, alphay, alphaz);
                amrex::ParallelFor(tex, tey, tez,
                [=] AMREX_GPU_DEVICE (int j, int k, int l)
                {
                    warpx_push_ex_f_ckc(j,k,l,Exfab,Ffab,
                                        betaxy, betaxz, betayx, betayz, betazx, betazy,
                                        gammax, gammay, gammaz,
                                        alphax, alphay, alphaz);
                },
                [=] AMREX_GPU_DEVICE (int j, int k, int l)
                {
                    warpx_push_ey_f_ckc(j,k,l,Eyfab,Ffab,
                                        betaxy, betaxz, betayx, betayz, betazx, betazy,
                                        gammax, gammay, gammaz,
                                        alphax, alphay, alphaz);
                },
                [=] AMREX_GPU_DEVICE (int j, int k, int l)
                {
                    warpx_push_ez_f_ckc(j,k,l,Ezfab,Ffab,
                                        betaxy, betaxz, betayx, betayz, betazx, betazy,
                                        gammax, gammay, gammaz,
                                        alphax, alphay, alphaz);
                });
            }
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
        if (F) pml[lev]->ExchangeF(patch_type, F, do_pml_in_domain);

        const auto& pml_B = (patch_type == PatchType::fine) ? pml[lev]->GetB_fp() : pml[lev]->GetB_cp();
        const auto& pml_E = (patch_type == PatchType::fine) ? pml[lev]->GetE_fp() : pml[lev]->GetE_cp();
        const auto& pml_j = (patch_type == PatchType::fine) ? pml[lev]->Getj_fp() : pml[lev]->Getj_cp();
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

            auto const& pml_Exfab = pml_E[0]->array(mfi);
            auto const& pml_Eyfab = pml_E[1]->array(mfi);
            auto const& pml_Ezfab = pml_E[2]->array(mfi);

            WRPX_PUSH_PML_EVEC(
			     tex.loVect(), tex.hiVect(),
			     tey.loVect(), tey.hiVect(),
			     tez.loVect(), tez.hiVect(),
			     BL_TO_FORTRAN_3D((*pml_E[0])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_E[1])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_E[2])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_B[0])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_B[1])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_B[2])[mfi]),
           &dtsdx_c2, &dtsdy_c2, &dtsdz_c2);

            if (pml_has_particles) {
                // Update the E field in the PML, using the current
                // deposited by the particles in the PML
                auto const& pml_jxfab = pml_j[0]->array(mfi);
                auto const& pml_jyfab = pml_j[1]->array(mfi);
                auto const& pml_jzfab = pml_j[2]->array(mfi);
                const Real* sigmaj_x = sigba[mfi].sigma[0].data();
                const Real* sigmaj_y = sigba[mfi].sigma[1].data();
                const Real* sigmaj_z = sigba[mfi].sigma[2].data();

                auto const& AMREX_RESTRICT x_lo = sigba[mfi].sigma[0].lo();
#if (AMREX_SPACEDIM == 3)
                auto const& AMREX_RESTRICT y_lo = sigba[mfi].sigma[1].lo();
                auto const& AMREX_RESTRICT z_lo = sigba[mfi].sigma[2].lo();
#else
                int y_lo = 0;
                auto const& AMREX_RESTRICT z_lo = sigba[mfi].sigma[1].lo();
#endif

                amrex::ParallelFor( tex, tey, tez,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                        push_ex_pml_current(i,j,k,
                            pml_Exfab, pml_jxfab, sigmaj_y, sigmaj_z,
                            y_lo, z_lo, mu_c2_dt);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                        push_ey_pml_current(i,j,k,
                            pml_Eyfab, pml_jyfab, sigmaj_x, sigmaj_z,
                            x_lo, z_lo, mu_c2_dt);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                        push_ez_pml_current(i,j,k,
                            pml_Ezfab, pml_jzfab, sigmaj_x, sigmaj_y,
                            x_lo, y_lo, mu_c2_dt);
                    }
                );
            }

            if (pml_F)
            {
                WRPX_PUSH_PML_EVEC_F(
				   tex.loVect(), tex.hiVect(),
				   tey.loVect(), tey.hiVect(),
				   tez.loVect(), tez.hiVect(),
				   BL_TO_FORTRAN_3D((*pml_E[0])[mfi]),
				   BL_TO_FORTRAN_3D((*pml_E[1])[mfi]),
				   BL_TO_FORTRAN_3D((*pml_E[2])[mfi]),
				   BL_TO_FORTRAN_3D((*pml_F   )[mfi]),
           &dtsdx_c2, &dtsdy_c2, &dtsdz_c2,
				   &WarpX::maxwell_fdtd_solver_id);
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

    const int patch_level = (patch_type == PatchType::fine) ? lev : lev-1;
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
            WRPX_PUSH_PML_F(bx.loVect(), bx.hiVect(),
			  BL_TO_FORTRAN_ANYD((*pml_F   )[mfi]),
			  BL_TO_FORTRAN_ANYD((*pml_E[0])[mfi]),
			  BL_TO_FORTRAN_ANYD((*pml_E[1])[mfi]),
			  BL_TO_FORTRAN_ANYD((*pml_E[2])[mfi]),
			  &dtsdx[0], &dtsdx[1], &dtsdx[2]);
        }
    }
}

#ifdef WARPX_DIM_RZ
// This scales the current by the inverse volume and wraps around the depostion at negative radius.
// It is faster to apply this on the grid than to do it particle by particle.
// It is put here since there isn't another nice place for it.
void
WarpX::ApplyInverseVolumeScalingToCurrentDensity (MultiFab* Jx, MultiFab* Jy, MultiFab* Jz, int lev)
{
    const long ngJ = Jx->nGrow();
    const std::array<Real,3>& dx = WarpX::CellSize(lev);
    const Real dr = dx[0];

    Box tilebox;

    for ( MFIter mfi(*Jx, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {

        Array4<Real> const& Jr_arr = Jx->array(mfi);
        Array4<Real> const& Jt_arr = Jy->array(mfi);
        Array4<Real> const& Jz_arr = Jz->array(mfi);

        tilebox = mfi.tilebox();
        Box tbr = convert(tilebox, WarpX::jx_nodal_flag);
        Box tbt = convert(tilebox, WarpX::jy_nodal_flag);
        Box tbz = convert(tilebox, WarpX::jz_nodal_flag);

        // Lower corner of tile box physical domain
        // Note that this is done before the tilebox.grow so that
        // these do not include the guard cells.
        const std::array<Real, 3>& xyzmin = WarpX::LowerCorner(tilebox, lev);
        const Dim3 lo = lbound(tilebox);
        const Real rmin = xyzmin[0];
        const int irmin = lo.x;

        // Rescale current in r-z mode since the inverse volume factor was not
        // included in the current deposition.
        amrex::ParallelFor(tbr, tbt, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // Wrap the current density deposited in the guard cells around
            // to the cells above the axis.
            // Note that Jr(i==0) is at 1/2 dr.
            if (rmin == 0. && 0 <= i && i < ngJ) {
                Jr_arr(i,j,0) -= Jr_arr(-1-i,j,0);
            }
            // Apply the inverse volume scaling
            // Since Jr is not cell centered in r, no need for distinction
            // between on axis and off-axis factors
            const amrex::Real r = std::abs(rmin + (i - irmin + 0.5)*dr);
            Jr_arr(i,j,0) /= (2.*MathConst::pi*r);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // Wrap the current density deposited in the guard cells around
            // to the cells above the axis.
            // Jt is located on the boundary
            if (rmin == 0. && 0 < i && i <= ngJ) {
                Jt_arr(i,j,0) += Jt_arr(-i,j,0);
            }

            // Apply the inverse volume scaling
            // Jt is forced to zero on axis.
            const amrex::Real r = std::abs(rmin + (i - irmin)*dr);
            if (r == 0.) {
                Jt_arr(i,j,0) = 0.;
            } else {
                Jt_arr(i,j,0) /= (2.*MathConst::pi*r);
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // Wrap the current density deposited in the guard cells around
            // to the cells above the axis.
            // Jz is located on the boundary
            if (rmin == 0. && 0 < i && i <= ngJ) {
                Jz_arr(i,j,0) += Jz_arr(-i,j,0);
            }

            // Apply the inverse volume scaling
            const amrex::Real r = std::abs(rmin + (i - irmin)*dr);
            if (r == 0.) {
                // Verboncoeur JCP 164, 421-427 (2001) : corrected volume on axis
                Jz_arr(i,j,0) /= (MathConst::pi*dr/3.);
            } else {
                Jz_arr(i,j,0) /= (2.*MathConst::pi*r);
            }
        });
    }
}

void
WarpX::ApplyInverseVolumeScalingToChargeDensity (MultiFab* Rho, int lev)
{
    const long ngRho = Rho->nGrow();
    const std::array<Real,3>& dx = WarpX::CellSize(lev);
    const Real dr = dx[0];

    Box tilebox;

    for ( MFIter mfi(*Rho, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {

        Array4<Real> const& Rho_arr = Rho->array(mfi);

        tilebox = mfi.tilebox();
        Box tb = convert(tilebox, IntVect::TheUnitVector());

        // Lower corner of tile box physical domain
        // Note that this is done before the tilebox.grow so that
        // these do not include the guard cells.
        const std::array<Real, 3>& xyzmin = WarpX::LowerCorner(tilebox, lev);
        const Dim3 lo = lbound(tilebox);
        const Real rmin = xyzmin[0];
        const int irmin = lo.x;

        // Rescale charge in r-z mode since the inverse volume factor was not
        // included in the charge deposition.
        amrex::ParallelFor(tb, Rho->nComp(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int icomp)
        {
            // Wrap the charge density deposited in the guard cells around
            // to the cells above the axis.
            // Rho is located on the boundary
            if (rmin == 0. && 0 < i && i <= ngRho) {
                Rho_arr(i,j,0,icomp) += Rho_arr(-i,j,0,icomp);
            }

            // Apply the inverse volume scaling
            const amrex::Real r = std::abs(rmin + (i - irmin)*dr);
            if (r == 0.) {
                // Verboncoeur JCP 164, 421-427 (2001) : corrected volume on axis
                Rho_arr(i,j,0,icomp) /= (MathConst::pi*dr/3.);
            } else {
                Rho_arr(i,j,0,icomp) /= (2.*MathConst::pi*r);
            }
        });
    }
}
#endif
