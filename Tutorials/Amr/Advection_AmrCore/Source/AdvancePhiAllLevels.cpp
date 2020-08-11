#include <AmrCoreAdv.H>
#include <Kernels.H>

#include <AMReX_MultiFabUtil.H>

using namespace amrex;

// advance all levels for a single time step
void
AmrCoreAdv::AdvancePhiAllLevels (Real time, Real dt_lev, int /*iteration*/)
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

        const auto dx = geom[lev].CellSizeArray();
        AMREX_D_TERM(Real dtdx = dt_lev/dx[0];,
                     Real dtdy = dt_lev/dx[1];,
                     Real dtdz = dt_lev/dx[2]);

        // State with ghost cells
        MultiFab Sborder(grids[lev], dmap[lev], phi_new[lev].nComp(), num_grow);
        FillPatch(lev, time, Sborder, 0, Sborder.nComp());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            FArrayBox tmpfab;
            for (MFIter mfi(phi_new[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                AMREX_D_TERM(Array4<Real const> velx = facevel[lev][0].const_array(mfi);,
                             Array4<Real const> vely = facevel[lev][1].const_array(mfi);,
                             Array4<Real const> velz = facevel[lev][2].const_array(mfi));

                const Box& bx = mfi.tilebox();
                const Box& gbx = amrex::grow(bx, 1);

                Array4<Real const> statein = Sborder.const_array(mfi);

                AMREX_D_TERM(Array4<Real> fluxx = fluxes[lev][0].array(mfi);,
                             Array4<Real> fluxy = fluxes[lev][1].array(mfi);,
                             Array4<Real> fluxz = fluxes[lev][2].array(mfi));

                int ntmpcomps = (AMREX_SPACEDIM == 2) ? 4 : 11;
                tmpfab.resize(amrex::grow(bx,2),ntmpcomps);
                Elixir tmpeli = tmpfab.elixir();
                int itmp = 0;

                Array4<Real> slope2 = tmpfab.array(itmp);
                Array4<Real const> slope2_c = slope2;
                itmp += 1;
                Array4<Real> slope4 = tmpfab.array(itmp);
                Array4<Real const> slope4_c = slope4;
                itmp += 1;

                // compute longitudinal fluxes
                // ===========================

                // x -------------------------
                Array4<Real> phix = tmpfab.array(itmp);
                Array4<Real const> phix_c = phix;
                itmp += 1;

                amrex::launch(amrex::grow(gbx,Direction::x,1),
                [=] AMREX_GPU_DEVICE (const Box& tbx)
                {
                    slopex2(tbx, slope2, statein);
                });

                amrex::launch(gbx,
                [=] AMREX_GPU_DEVICE (const Box& tbx)
                {
                    slopex4(tbx, slope4, statein, slope2_c);
                });

                Box b = gbx;
                amrex::ParallelFor(b.grow(Direction::x,-1).surroundingNodes(Direction::x),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_x(i, j, k, phix, statein, velx, slope4_c, dtdx);
                });

                // y -------------------------
                Array4<Real> phiy = tmpfab.array(itmp);
                Array4<Real const> phiy_c = phiy;
                itmp += 1;

                amrex::launch(amrex::grow(gbx,Direction::y,1),
                [=] AMREX_GPU_DEVICE (const Box& tbx)
                {
                    slopey2(tbx, slope2, statein);
                });

                amrex::launch(gbx,
                [=] AMREX_GPU_DEVICE (const Box& tbx)
                {
                    slopey4(tbx, slope4, statein, slope2_c);
                });

                b = gbx;
                amrex::ParallelFor(b.grow(Direction::y,-1).surroundingNodes(Direction::y),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_y(i, j, k, phiy, statein, vely, slope4_c, dtdy);
                });

#if (AMREX_SPACEDIM > 2)
                // z -------------------------
                Array4<Real> phiz = tmpfab.array(itmp);
                Array4<Real const> phiz_c = phiz;
                itmp += 1;

                amrex::launch(amrex::grow(gbx,Direction::z,1),
                [=] AMREX_GPU_DEVICE (const Box& tbx)
                {
                    slopez2(tbx, slope2, statein);
                });

                amrex::launch(gbx,
                [=] AMREX_GPU_DEVICE (const Box& tbx)
                {
                    slopez4(tbx, slope4, statein, slope2_c);
                });

                b = gbx;
                amrex::ParallelFor(b.grow(Direction::z,-1).surroundingNodes(Direction::z),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_z(i, j, k, phiz, statein, velz, slope4_c, dtdz);
                });

                // compute transverse fluxes (3D only)
                // ===================================

                // xy --------------------
                Array4<Real> phix_y = tmpfab.array(itmp);
                Array4<Real const> phix_y_c = phix_y;
                itmp += 1;

                b = bx;
                amrex::ParallelFor(b.grow(Direction::z,1).surroundingNodes(Direction::x),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_xy(i, j, k, phix_y,
                            velx, vely,
                            phix_c, phiy_c,
                            dtdy);
                }); 

                // xz --------------------
                Array4<Real> phix_z = tmpfab.array(itmp);
                Array4<Real const> phix_z_c = phix_z;
                itmp += 1;

                b = bx;
                amrex::ParallelFor(b.grow(Direction::y,1).surroundingNodes(Direction::x),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_xz(i, j, k, phix_z,
                            velx, velz,
                            phix, phiz,
                            dtdz);
                }); 

                // yx --------------------
                Array4<Real> phiy_x = tmpfab.array(itmp);
                Array4<Real const> phiy_x_c = phiy_x;
                itmp += 1;

                b = bx;
                amrex::ParallelFor(b.grow(Direction::z,1).surroundingNodes(Direction::y),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_yx(i, j, k, phiy_x,
                            velx, vely,
                            phix, phiy,
                            dtdx);
                }); 

                // yz --------------------
                Array4<Real> phiy_z = tmpfab.array(itmp);
                Array4<Real const> phiy_z_c = phiy_z;
                itmp += 1;

                b = bx;
                amrex::ParallelFor(b.grow(Direction::x,1).surroundingNodes(Direction::y),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_yz(i, j, k, phiy_z,
                            vely, velz,
                            phiy, phiz,
                            dtdz);
                }); 

                // zx --------------------
                Array4<Real> phiz_x = tmpfab.array(itmp);
                Array4<Real const> phiz_x_c = phiz_x;
                itmp += 1;

                b = bx;
                amrex::ParallelFor(b.grow(Direction::y,1).surroundingNodes(Direction::z),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_zx(i, j, k, phiz_x,
                            velx, velz,
                            phix, phiz,
                            dtdx);
                }); 

                // zy --------------------
                Array4<Real> phiz_y = tmpfab.array(itmp);
                Array4<Real const> phiz_y_c = phiz_y;
                itmp += 1;

                b = bx;
                amrex::ParallelFor(b.grow(Direction::x,1).surroundingNodes(Direction::z),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    flux_zy(i, j, k, phiz_y,
                            vely, velz,
                            phiy, phiz,
                            dtdy);
                }); 
#endif

                // final edge states 
                // ===========================
                amrex::ParallelFor(mfi.nodaltilebox(0),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    create_flux_x(i, j, k, fluxx,
                                  AMREX_D_DECL(velx,vely,velz),
#if (AMREX_SPACEDIM == 3)
                                  phix_c, phiy_z_c, phiz_y_c,
                                  dtdy, dtdz);
#else
                                  phix_c, phiy_c,
                                  dtdy);
#endif
                });

                amrex::ParallelFor(mfi.nodaltilebox(1),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    create_flux_y(i, j, k, fluxy,
                                  AMREX_D_DECL(velx,vely,velz),
#if (AMREX_SPACEDIM == 3)
                                  phiy_c, phix_z_c, phiz_x_c,
                                  dtdx, dtdz);
#else
                                  phiy_c, phix_c,
                                  dtdx);
#endif
                });

#if (AMREX_SPACEDIM == 3)
                amrex::ParallelFor(mfi.nodaltilebox(2),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    create_flux_z(i, j, k, fluxz,
                                  velx, vely, velz,
                                  phiz_c, phix_y_c, phiy_x_c,
                                  dtdx, dtdy);
                });
#endif
                AMREX_ASSERT(itmp == ntmpcomps);
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
            AMREX_D_TERM(Real dtdx = dt_lev/dx[0];,
                         Real dtdy = dt_lev/dx[1];,
                         Real dtdz = dt_lev/dx[2]);

            // ===========================================
            // Compute phi_new using a conservative update 
            // ===========================================
            for (MFIter mfi(phi_new[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Array4<Real const> statein  = phi_old[lev].const_array(mfi);
                Array4<Real      > stateout = phi_new[lev].array(mfi);

                AMREX_D_TERM(Array4<Real const> fluxx = fluxes[lev][0].const_array(mfi);,
                             Array4<Real const> fluxy = fluxes[lev][1].const_array(mfi);,
                             Array4<Real const> fluxz = fluxes[lev][2].const_array(mfi));

                const Box& bx  = mfi.tilebox();
    
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    conservative(i, j, k,
                                 stateout, statein,
                                 AMREX_D_DECL(fluxx,fluxy,fluxz),
                                 AMREX_D_DECL(dtdx,dtdy,dtdz));
                });
            } // end mfi
        } // end omp
    } // end lev
}
