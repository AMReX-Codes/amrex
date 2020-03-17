/* Copyright 2019 Andrew Myers, Aurore Blelly, Axel Huebl
 * David Grote, Maxence Thevenet, Remi Lehe
 * Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "Utils/WarpXConst.H"
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
WarpX::EvolveB (amrex::Real a_dt)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        EvolveB(lev, a_dt);
    }
}

void
WarpX::EvolveB (int lev, amrex::Real a_dt)
{
    WARPX_PROFILE("WarpX::EvolveB()");
    EvolveB(lev, PatchType::fine, a_dt);
    if (lev > 0)
    {
        EvolveB(lev, PatchType::coarse, a_dt);
    }
}

void
WarpX::EvolveB (int lev, PatchType patch_type, amrex::Real a_dt)
{

    if (patch_type == PatchType::fine) {
        m_fdtd_solver_fp[lev]->EvolveB( Bfield_fp[lev], Efield_fp[lev], a_dt );
    } else {
        m_fdtd_solver_cp[lev]->EvolveB( Bfield_cp[lev], Efield_cp[lev], a_dt );
    }

    const int patch_level = (patch_type == PatchType::fine) ? lev : lev-1;
    const std::array<Real,3>& dx = WarpX::CellSize(patch_level);
    const Real dtsdx = a_dt/dx[0], dtsdy = a_dt/dx[1], dtsdz = a_dt/dx[2];

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
            auto const& pml_Bxfab = pml_B[0]->array(mfi);
            auto const& pml_Byfab = pml_B[1]->array(mfi);
            auto const& pml_Bzfab = pml_B[2]->array(mfi);
            auto const& pml_Exfab = pml_E[0]->array(mfi);
            auto const& pml_Eyfab = pml_E[1]->array(mfi);
            auto const& pml_Ezfab = pml_E[2]->array(mfi);
            if (WarpX::maxwell_fdtd_solver_id == 0) {
               amrex::ParallelFor(tbx, tby, tbz,
               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                   warpx_push_pml_bx_yee(i,j,k,pml_Bxfab,pml_Eyfab,pml_Ezfab,
                                        dtsdy,dtsdz);
               },
               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                   warpx_push_pml_by_yee(i,j,k,pml_Byfab,pml_Exfab,pml_Ezfab,
                                         dtsdx,dtsdz);
               },
               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                   warpx_push_pml_bz_yee(i,j,k,pml_Bzfab,pml_Exfab,pml_Eyfab,
                                        dtsdx,dtsdy);
               });
            }  else if (WarpX::maxwell_fdtd_solver_id == 1) {
               Real betaxy, betaxz, betayx, betayz, betazx, betazy;
               Real gammax, gammay, gammaz;
               Real alphax, alphay, alphaz;
               warpx_calculate_ckc_coefficients(dtsdx, dtsdy, dtsdz,
                                                betaxy, betaxz, betayx, betayz,
                                                betazx, betazy, gammax, gammay,
                                                gammaz, alphax, alphay, alphaz);

               amrex::ParallelFor(tbx, tby, tbz,
               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                   warpx_push_pml_bx_ckc(i,j,k,pml_Bxfab,pml_Eyfab,pml_Ezfab,
                                         betaxy, betaxz, betayx, betayz,
                                         betazx, betazy, gammax, gammay,
                                         gammaz, alphax, alphay, alphaz);
               },
               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                   warpx_push_pml_by_ckc(i,j,k,pml_Byfab,pml_Exfab,pml_Ezfab,
                                         betaxy, betaxz, betayx, betayz,
                                         betazx, betazy, gammax, gammay,
                                         gammaz, alphax, alphay, alphaz);
               },
               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                   warpx_push_pml_bz_ckc(i,j,k,pml_Bzfab,pml_Exfab,pml_Eyfab,
                                         betaxy, betaxz, betayx, betayz,
                                         betazx, betazy, gammax, gammay,
                                         gammaz, alphax, alphay, alphaz);
               });

            }
        }
    }
}

void
WarpX::EvolveE (amrex::Real a_dt)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        EvolveE(lev, a_dt);
    }
}

void
WarpX::EvolveE (int lev, amrex::Real a_dt)
{
    WARPX_PROFILE("WarpX::EvolveE()");
    EvolveE(lev, PatchType::fine, a_dt);
    if (lev > 0)
    {
        EvolveE(lev, PatchType::coarse, a_dt);
    }
}

void
WarpX::EvolveE (int lev, PatchType patch_type, amrex::Real a_dt)
{

    if (patch_type == PatchType::fine) {
        m_fdtd_solver_fp[lev]->EvolveE( Efield_fp[lev], Bfield_fp[lev],
                                      current_fp[lev], F_fp[lev], a_dt );
    } else {
        m_fdtd_solver_cp[lev]->EvolveE( Efield_cp[lev], Bfield_cp[lev],
                                      current_cp[lev], F_cp[lev], a_dt );
    }

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

    MultiFab* cost = WarpX::getCosts(lev);
    const IntVect& rr = (lev > 0) ? refRatio(lev-1) : IntVect::TheUnitVector();

    // xmin is only used by the kernel for cylindrical geometry,
    // in which case it is actually rmin.
    const Real xmin = Geom(0).ProbLo(0);

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
            auto const& pml_Bxfab = pml_B[0]->array(mfi);
            auto const& pml_Byfab = pml_B[1]->array(mfi);
            auto const& pml_Bzfab = pml_B[2]->array(mfi);

            amrex::ParallelFor(tex, tey, tez,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                warpx_push_pml_ex_yee(i,j,k,pml_Exfab,pml_Byfab,pml_Bzfab,
                                      dtsdy_c2,dtsdz_c2);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                warpx_push_pml_ey_yee(i,j,k,pml_Eyfab,pml_Bxfab,pml_Bzfab,
                                      dtsdx_c2,dtsdz_c2);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                warpx_push_pml_ez_yee(i,j,k,pml_Ezfab,pml_Bxfab,pml_Byfab,
                                      dtsdx_c2,dtsdy_c2);
            });

            if (pml_has_particles) {
                // Update the E field in the PML, using the current
                // deposited by the particles in the PML
                auto const& pml_jxfab = pml_j[0]->array(mfi);
                auto const& pml_jyfab = pml_j[1]->array(mfi);
                auto const& pml_jzfab = pml_j[2]->array(mfi);
                const Real* sigmaj_x = sigba[mfi].sigma[0].data();
                const Real* sigmaj_y = sigba[mfi].sigma[1].data();
                const Real* sigmaj_z = sigba[mfi].sigma[2].data();

                int const x_lo = sigba[mfi].sigma[0].lo();
#if (AMREX_SPACEDIM == 3)
                int const y_lo = sigba[mfi].sigma[1].lo();
                int const z_lo = sigba[mfi].sigma[2].lo();
#else
                int const y_lo = 0;
                int const z_lo = sigba[mfi].sigma[1].lo();
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

               auto const& pml_F_fab = pml_F->array(mfi);

               if (WarpX::maxwell_fdtd_solver_id == 0) {

                  amrex::ParallelFor(tex, tey, tez,
                  [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                      warpx_push_pml_ex_f_yee(i,j,k,pml_Exfab,pml_F_fab,dtsdx_c2);
                  },
                  [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                      warpx_push_pml_ey_f_yee(i,j,k,pml_Eyfab,pml_F_fab,dtsdy_c2);
                  },
                  [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                      warpx_push_pml_ez_f_yee(i,j,k,pml_Ezfab,pml_F_fab,dtsdz_c2);
                  });

               } else if (WarpX::maxwell_fdtd_solver_id == 1) {

                  Real betaxy, betaxz, betayx, betayz, betazx, betazy;
                  Real gammax, gammay, gammaz;
                  Real alphax, alphay, alphaz;
                  warpx_calculate_ckc_coefficients(dtsdx_c2, dtsdy_c2, dtsdz_c2,
                                                   betaxy, betaxz, betayx, betayz,
                                                   betazx, betazy, gammax, gammay,
                                                   gammaz, alphax, alphay, alphaz);
                  amrex::ParallelFor(tex, tey, tez,
                  [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                      warpx_push_pml_ex_f_ckc(i,j,k,pml_Exfab,pml_F_fab,
                                              alphax,betaxy,betaxz,gammax);
                  },
                  [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                      warpx_push_pml_ey_f_ckc(i,j,k,pml_Eyfab,pml_F_fab,
                                              alphay,betayx,betayz,gammay);
                  },
                  [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                      warpx_push_pml_ez_f_ckc(i,j,k,pml_Ezfab,pml_F_fab,
                                              alphaz,betazx,betazy,gammaz);
                  });

               }
            }
        }
    }
}

void
WarpX::EvolveF (amrex::Real a_dt, DtType a_dt_type)
{
    if (!do_dive_cleaning) return;

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        EvolveF(lev, a_dt, a_dt_type);
    }
}

void
WarpX::EvolveF (int lev, amrex::Real a_dt, DtType a_dt_type)
{
    if (!do_dive_cleaning) return;

    EvolveF(lev, PatchType::fine, a_dt, a_dt_type);
    if (lev > 0) EvolveF(lev, PatchType::coarse, a_dt, a_dt_type);
}

void
WarpX::EvolveF (int lev, PatchType patch_type, amrex::Real a_dt, DtType a_dt_type)
{
    if (!do_dive_cleaning) return;

    WARPX_PROFILE("WarpX::EvolveF()");

    const int rhocomp = (a_dt_type == DtType::FirstHalf) ? 0 : 1;

    if (patch_type == PatchType::fine) {
        m_fdtd_solver_fp[lev]->EvolveF( F_fp[lev], Efield_fp[lev],
                                        rho_fp[lev], rhocomp, a_dt );
    } else {
        m_fdtd_solver_cp[lev]->EvolveF( F_cp[lev], Efield_cp[lev],
                                        rho_cp[lev], rhocomp, a_dt );
    }

    const int patch_level = (patch_type == PatchType::fine) ? lev : lev-1;
    const auto& dx = WarpX::CellSize(patch_level);
    const std::array<Real,3> dtsdx {a_dt/dx[0], a_dt/dx[1], a_dt/dx[2]};

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

    constexpr int NODE = amrex::IndexType::NODE;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Jx->ixType().ixType()[0] != NODE,
        "Jr should never node-centered in r");

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
        std::array<amrex::Real,3> galilean_shift = {0,0,0};
        const std::array<Real, 3>& xyzmin = WarpX::LowerCorner(tilebox, galilean_shift, lev);
        const Real rmin  = xyzmin[0] + (tbr.type(0) == NODE ? 0. : 0.5*dx[0]);
        const Real rminr = xyzmin[0] + (tbr.type(0) == NODE ? 0. : 0.5*dx[0]);
        const Real rmint = xyzmin[0] + (tbt.type(0) == NODE ? 0. : 0.5*dx[0]);
        const Real rminz = xyzmin[0] + (tbz.type(0) == NODE ? 0. : 0.5*dx[0]);
        const Dim3 lo = lbound(tilebox);
        const int irmin = lo.x;
        int const ishift_t = (rmint > rmin ? 1 : 0);
        int const ishift_z = (rminz > rmin ? 1 : 0);

        const long nmodes = n_rz_azimuthal_modes;

        // Grow the tileboxes to include the guard cells, except for the
        // guard cells at negative radius.
        if (rmin > 0.) {
           tbr.growLo(0, ngJ);
           tbt.growLo(0, ngJ);
           tbz.growLo(0, ngJ);
        }
        tbr.growHi(0, ngJ);
        tbt.growHi(0, ngJ);
        tbz.growHi(0, ngJ);
        tbr.grow(1, ngJ);
        tbt.grow(1, ngJ);
        tbz.grow(1, ngJ);

        // Rescale current in r-z mode since the inverse volume factor was not
        // included in the current deposition.
        amrex::ParallelFor(tbr, tbt, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // Wrap the current density deposited in the guard cells around
            // to the cells above the axis.
            // Note that Jr(i==0) is at 1/2 dr.
            if (rmin == 0. && 0 <= i && i < ngJ) {
                Jr_arr(i,j,0,0) -= Jr_arr(-1-i,j,0,0);
            }
            // Apply the inverse volume scaling
            // Since Jr is never node centered in r, no need for distinction
            // between on axis and off-axis factors
            const amrex::Real r = std::abs(rminr + (i - irmin)*dr);
            Jr_arr(i,j,0,0) /= (2.*MathConst::pi*r);

            for (int imode=1 ; imode < nmodes ; imode++) {
                const Real ifact = ( (imode%2) == 0 ? +1. : -1.);
                // Wrap the current density deposited in the guard cells around
                // to the cells above the axis.
                // Note that Jr(i==0) is at 1/2 dr.
                if (rmin == 0. && 0 <= i && i < ngJ) {
                    Jr_arr(i,j,0,2*imode-1) -= ifact*Jr_arr(-1-i,j,0,2*imode-1);
                    Jr_arr(i,j,0,2*imode) -= ifact*Jr_arr(-1-i,j,0,2*imode);
                }
                // Apply the inverse volume scaling
                // Since Jr is never node centered in r, no need for distinction
                // between on axis and off-axis factors
                Jr_arr(i,j,0,2*imode-1) /= (2.*MathConst::pi*r);
                Jr_arr(i,j,0,2*imode) /= (2.*MathConst::pi*r);
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // Wrap the current density deposited in the guard cells around
            // to the cells above the axis.
            // If Jt is node centered, Jt[0] is located on the boundary.
            // If Jt is cell centered, Jt[0] is at 1/2 dr.
            if (rmin == 0. && 0 < i && i <= ngJ-ishift_t) {
                Jt_arr(i,j,0,0) += Jt_arr(-ishift_t-i,j,0,0);
            }

            // Apply the inverse volume scaling
            // Jt is forced to zero on axis.
            const amrex::Real r = std::abs(rmint + (i - irmin)*dr);
            if (r == 0.) {
                Jt_arr(i,j,0,0) = 0.;
            } else {
                Jt_arr(i,j,0,0) /= (2.*MathConst::pi*r);
            }

            for (int imode=1 ; imode < nmodes ; imode++) {
                const Real ifact = ( (imode%2) == 0 ? +1. : -1.);
                // Wrap the current density deposited in the guard cells around
                // to the cells above the axis.
                if (rmin == 0. && 0 < i && i <= ngJ-ishift_t) {
                    Jt_arr(i,j,0,2*imode-1) += ifact*Jt_arr(-ishift_t-i,j,0,2*imode-1);
                    Jt_arr(i,j,0,2*imode) += ifact*Jt_arr(-ishift_t-i,j,0,2*imode);
                }

                // Apply the inverse volume scaling
                // Jt is forced to zero on axis.
                if (r == 0.) {
                    Jt_arr(i,j,0,2*imode-1) = 0.;
                    Jt_arr(i,j,0,2*imode) = 0.;
                } else {
                    Jt_arr(i,j,0,2*imode-1) /= (2.*MathConst::pi*r);
                    Jt_arr(i,j,0,2*imode) /= (2.*MathConst::pi*r);
                }
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // Wrap the current density deposited in the guard cells around
            // to the cells above the axis.
            // If Jz is node centered, Jt[0] is located on the boundary.
            // If Jz is cell centered, Jt[0] is at 1/2 dr.
            if (rmin == 0. && 0 < i && i <= ngJ-ishift_z) {
                Jz_arr(i,j,0,0) += Jz_arr(-ishift_z-i,j,0,0);
            }

            // Apply the inverse volume scaling
            const amrex::Real r = std::abs(rminz + (i - irmin)*dr);
            if (r == 0.) {
                // Verboncoeur JCP 164, 421-427 (2001) : corrected volume on axis
                Jz_arr(i,j,0,0) /= (MathConst::pi*dr/3.);
            } else {
                Jz_arr(i,j,0,0) /= (2.*MathConst::pi*r);
            }

            for (int imode=1 ; imode < nmodes ; imode++) {
                const Real ifact = ( (imode%2) == 0 ? +1. : -1.);
                // Wrap the current density deposited in the guard cells around
                // to the cells above the axis.
                if (rmin == 0. && 0 < i && i <= ngJ-ishift_z) {
                    Jz_arr(i,j,0,2*imode-1) += ifact*Jz_arr(-ishift_z-i,j,0,2*imode-1);
                    Jz_arr(i,j,0,2*imode) += ifact*Jz_arr(-ishift_z-i,j,0,2*imode);
                }

                // Apply the inverse volume scaling
                if (r == 0.) {
                    // Verboncoeur JCP 164, 421-427 (2001) : corrected volume on axis
                    Jz_arr(i,j,0,2*imode-1) /= (MathConst::pi*dr/3.);
                    Jz_arr(i,j,0,2*imode) /= (MathConst::pi*dr/3.);
                } else {
                    Jz_arr(i,j,0,2*imode-1) /= (2.*MathConst::pi*r);
                    Jz_arr(i,j,0,2*imode) /= (2.*MathConst::pi*r);
                }
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

    constexpr int NODE = amrex::IndexType::NODE;

    Box tilebox;

    for ( MFIter mfi(*Rho, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {

        Array4<Real> const& Rho_arr = Rho->array(mfi);

        tilebox = mfi.tilebox();
        Box tb = convert(tilebox, rho_nodal_flag);

        // Lower corner of tile box physical domain
        // Note that this is done before the tilebox.grow so that
        // these do not include the guard cells.
        std::array<amrex::Real,3> galilean_shift = {0,0,0};
        const std::array<Real, 3>& xyzmin = WarpX::LowerCorner(tilebox, galilean_shift, lev);
        const Dim3 lo = lbound(tilebox);
        const Real rmin = xyzmin[0];
        const Real rminr = xyzmin[0] + (tb.type(0) == NODE ? 0. : 0.5*dx[0]);
        const int irmin = lo.x;
        int ishift = (rminr > rmin ? 1 : 0);

        // Grow the tilebox to include the guard cells, except for the
        // guard cells at negative radius.
        if (rmin > 0.) {
           tb.growLo(0, ngRho);
        }
        tb.growHi(0, ngRho);
        tb.grow(1, ngRho);

        // Rescale charge in r-z mode since the inverse volume factor was not
        // included in the charge deposition.
        // Note that the loop is also over ncomps, which takes care of the RZ modes,
        // as well as the old and new rho.
        amrex::ParallelFor(tb, Rho->nComp(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int icomp)
        {
            // Wrap the charge density deposited in the guard cells around
            // to the cells above the axis.
            // Rho is located on the boundary
            if (rmin == 0. && 0 < i && i <= ngRho-ishift) {
                Rho_arr(i,j,0,icomp) += Rho_arr(-ishift-i,j,0,icomp);
            }

            // Apply the inverse volume scaling
            const amrex::Real r = std::abs(rminr + (i - irmin)*dr);
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
