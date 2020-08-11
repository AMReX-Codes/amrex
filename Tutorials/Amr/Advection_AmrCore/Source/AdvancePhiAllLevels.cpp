#include <AmrCoreAdv.H>
#include <Kernels.H>

#include <AMReX_MultiFabUtil.H>

using namespace amrex;

// advance all levels for a single time step
void
AmrCoreAdv::AdvancePhiAllLevels (Real time, Real dt_lev, int iteration)
{
    constexpr int num_grow = 3;

    Vector< Array<MultiFab,AMREX_SPACEDIM> > fluxes(finest_level+1);
    for (int lev = 0; lev <= finest_level; lev++)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            BoxArray ba = grids[lev];
            ba.surroundingNodes(idim);
            fluxes[lev][idim] = MultiFab(ba, dmap[lev], 1, 0);
        }
    }

    for (int lev = 0; lev <= finest_level; lev++)
    {
        std::swap(phi_old[lev], phi_new[lev]);
        t_old[lev] = t_new[lev];
        t_new[lev] += dt_lev;

        const Real old_time = t_old[lev];
        const Real new_time = t_new[lev];
        const Real ctr_time = 0.5*(old_time+new_time);

        const auto dx = geom[lev].CellSizeArray();
        GpuArray<Real, AMREX_SPACEDIM> dtdx;
        for (int i=0; i<AMREX_SPACEDIM; ++i)
            dtdx[i] = dt_lev/(dx[i]);

        const Real* prob_lo = geom[lev].ProbLo();

        // State with ghost cells
        MultiFab Sborder(grids[lev], dmap[lev], phi_new[lev].nComp(), num_grow);
        FillPatch(lev, time, Sborder, 0, Sborder.nComp());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            for (MFIter mfi(phi_new[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                GpuArray<Array4<Real>, AMREX_SPACEDIM> vel{ AMREX_D_DECL( facevel[lev][0].array(mfi),
                                                                          facevel[lev][1].array(mfi),
                                                                          facevel[lev][2].array(mfi)) };

                const Box& bx = mfi.tilebox();
                const Box& gbx = amrex::grow(bx, 1);

                Array4<Real> statein  = Sborder.array(mfi);
                Array4<Real> stateout = phi_new[lev].array(mfi);

                GpuArray<Array4<Real>, AMREX_SPACEDIM> flux{ AMREX_D_DECL(fluxes[lev][0].array(mfi),
                                                                          fluxes[lev][1].array(mfi),
                                                                          fluxes[lev][2].array(mfi)) };
    
                    AMREX_D_TERM(const Box& dqbxx = amrex::grow(bx, IntVect{AMREX_D_DECL(2, 1, 1)});,
                             const Box& dqbxy = amrex::grow(bx, IntVect{AMREX_D_DECL(1, 2, 1)});,
                             const Box& dqbxz = amrex::grow(bx, IntVect{AMREX_D_DECL(1, 1, 2)}););

                FArrayBox slope2fab (amrex::grow(bx, 2), 1);
                Elixir slope2eli = slope2fab.elixir();
                Array4<Real> slope2 = slope2fab.array();
                FArrayBox slope4fab (amrex::grow(bx, 1), 1);
                Elixir slope4eli = slope4fab.elixir();
                Array4<Real> slope4 = slope4fab.array();
    
                // compute longitudinal fluxes
                // ===========================

                // x -------------------------
                FArrayBox phixfab (gbx, 1);
                Elixir phixeli = phixfab.elixir();
                Array4<Real> phix = phixfab.array();

                amrex::launch(dqbxx,
                [=] AMREX_GPU_DEVICE (const Box& tbx)
                {
                    slopex2(tbx, statein, slope2);
                });

                amrex::launch(gbx,
                [=] AMREX_GPU_DEVICE (const Box& tbx)
                {
                    slopex4(tbx, statein, slope2, slope4);
                });

                amrex::ParallelFor(amrex::growLo(gbx, 0, -1),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_x(i, j, k, statein, vel[0], phix, slope4, dtdx); 
                });

                // y -------------------------
                FArrayBox phiyfab (gbx, 1);
                Elixir phiyeli = phiyfab.elixir();
                Array4<Real> phiy = phiyfab.array();

                amrex::launch(dqbxy,
                [=] AMREX_GPU_DEVICE (const Box& tbx)
                {
                    slopey2(tbx, statein, slope2);
                });

                amrex::launch(gbx,
                [=] AMREX_GPU_DEVICE (const Box& tbx)
                {
                    slopey4(tbx, statein, slope2, slope4);
                });

                amrex::ParallelFor(amrex::growLo(gbx, 1, -1),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_y(i, j, k, statein, vel[1], phiy, slope4, dtdx); 
                });

#if (AMREX_SPACEDIM > 2)
                // z -------------------------
                FArrayBox phizfab (gbx, 1);
                Elixir phizeli = phizfab.elixir();
                Array4<Real> phiz = phizfab.array();

                amrex::launch(dqbxz,
                [=] AMREX_GPU_DEVICE (const Box& tbx)
                {
                    slopez2(tbx, statein, slope2);
                });

                amrex::launch(gbx,
                [=] AMREX_GPU_DEVICE (const Box& tbx)
                {
                    slopez4(tbx, statein, slope2, slope4);
                });

                amrex::ParallelFor(amrex::growLo(gbx, 2, -1),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_z(i, j, k, statein, vel[2], phiz, slope4, dtdx); 
                });

                // compute transverse fluxes (3D only)
                // ===================================

                AMREX_D_TERM(const Box& gbxx = amrex::grow(bx, 0, 1);,
                             const Box& gbxy = amrex::grow(bx, 1, 1);,
                             const Box& gbxz = amrex::grow(bx, 2, 1););

                // xy --------------------
                FArrayBox phix_yfab (gbx, 1);
                Elixir phix_yeli = phix_yfab.elixir();
                Array4<Real> phix_y = phix_yfab.array();
    
                amrex::ParallelFor(amrex::growHi(gbxz, 0, 1),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_xy(i, j, k, 
                            AMREX_D_DECL(vel[0], vel[1], vel[2]),
                            AMREX_D_DECL(phix, phiy, phiz),
                            phix_y, dtdx);
                }); 

                // xz --------------------
                FArrayBox phix_zfab (gbx, 1);
                Elixir phix_zeli = phix_zfab.elixir();
                Array4<Real> phix_z = phix_zfab.array();

                amrex::ParallelFor(amrex::growHi(gbxy, 0, 1),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_xz(i, j, k,
                            AMREX_D_DECL(vel[0], vel[1], vel[2]),
                            AMREX_D_DECL(phix, phiy, phiz),
                            phix_z, dtdx);
                }); 

                // yx --------------------
                FArrayBox phiy_xfab (gbx, 1);
                FArrayBox phiy_zfab (gbx, 1);
                Elixir phiy_xeli = phiy_xfab.elixir();
                Elixir phiy_zeli = phiy_zfab.elixir();
                Array4<Real> phiy_x = phiy_xfab.array();
                Array4<Real> phiy_z = phiy_zfab.array();

                amrex::ParallelFor(amrex::growHi(gbxz, 1, 1),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_yx(i, j, k,
                            AMREX_D_DECL(vel[0], vel[1], vel[2]),
                            AMREX_D_DECL(phix, phiy, phiz),
                            phiy_x, dtdx);
                }); 

                // yz --------------------
                amrex::ParallelFor(amrex::growHi(gbxx, 1, 1),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_yz(i, j, k,
                            AMREX_D_DECL(vel[0], vel[1], vel[2]),
                            AMREX_D_DECL(phix, phiy, phiz),
                            phiy_z, dtdx);
                }); 

                // zx & zy --------------------
                FArrayBox phiz_xfab (gbx, 1);
                FArrayBox phiz_yfab (gbx, 1);
                Elixir phiz_xeli = phiz_xfab.elixir();
                Elixir phiz_yeli = phiz_yfab.elixir();
                Array4<Real> phiz_x = phiz_xfab.array();
                Array4<Real> phiz_y = phiz_yfab.array();
    
                amrex::ParallelFor(amrex::growHi(gbxy, 2, 1),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_zx(i, j, k, 
                            AMREX_D_DECL(vel[0], vel[1], vel[2]),
                            AMREX_D_DECL(phix, phiy, phiz),
                            phiz_x, dtdx);
                }); 

                amrex::ParallelFor(amrex::growHi(gbxx, 2, 1),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_zy(i, j, k,
                            AMREX_D_DECL(vel[0], vel[1], vel[2]),
                            AMREX_D_DECL(phix, phiy, phiz),
                            phiz_y, dtdx);
                }); 
#endif

                // final edge states 
                // ===========================
                amrex::ParallelFor(amrex::growHi(bx, 0, 1),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    create_flux_x(i, j, k,
                                  vel[0], vel[1], 
#if (AMREX_SPACEDIM > 2)
                                  vel[2],
#endif
#if (AMREX_SPACEDIM > 2)
                                  phix, phiy_z, phiz_y,
#else
                                  phix, phiy,
#endif
                                  flux[0], dtdx);
                });

                amrex::ParallelFor(amrex::growHi(bx, 1, 1),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    create_flux_y(i, j, k,
                                  vel[0], vel[1], 
#if (AMREX_SPACEDIM > 2)
                                  vel[2],
#endif
#if (AMREX_SPACEDIM > 2)
                                  phiy, phix_z, phiz_x,
#else
                                  phiy, phix,
#endif
                                  flux[1], dtdx);
            });

#if (AMREX_SPACEDIM > 2)
            amrex::ParallelFor(amrex::growHi(bx, 2, 1),
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                create_flux_z(i, j, k,
                               vel[0], vel[1], vel[2],
                               phiz, phix_y, phiy_x,
                               flux[2], dtdx);
            });
#endif
            } // end mfi
        } // end omp
    } // end lev

    // =======================================================
    // Average down the fluxes before using them to update phi 
    // =======================================================
    for (int lev = finest_level; lev > 0; lev--)
    {
       average_down_faces(amrex::GetArrOfConstPtrs(fluxes[lev  ]),
                          amrex::GetArrOfPtrs     (fluxes[lev-1]),
                          refRatio(lev-1), Geom(lev-1));
    } 

    for (int lev = 0; lev <= finest_level; lev++)
    {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            const auto dx = geom[lev].CellSizeArray();
            GpuArray<Real, AMREX_SPACEDIM> dtdx;
            for (int i=0; i<AMREX_SPACEDIM; ++i)
                dtdx[i] = dt_lev/(dx[i]);

            // ===========================================
            // Compute phi_new using a conservative update 
            // ===========================================
            for (MFIter mfi(phi_new[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Array4<Real> statein  = phi_old[lev].array(mfi);
                Array4<Real> stateout = phi_new[lev].array(mfi);

                GpuArray<Array4<Real>, AMREX_SPACEDIM> flux{ AMREX_D_DECL(fluxes[lev][0].array(mfi),
                                                                          fluxes[lev][1].array(mfi),
                                                                          fluxes[lev][2].array(mfi)) };

                const Box& bx  = mfi.tilebox();
    
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    conservative(i, j, k,
                                 statein, stateout,
                                 AMREX_D_DECL(flux[0], flux[1], flux[2]),
                                 dtdx);
                });
            } // end mfi
        } // end omp
    } // end lev
}
