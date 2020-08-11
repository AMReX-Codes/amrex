#include <AmrCoreAdv.H>
#include <Kernels.H>

using namespace amrex;

// Advance a single level for a single time step, updates flux registers
void
AmrCoreAdv::AdvancePhiAtLevel (int lev, Real time, Real dt_lev, int iteration, int ncycle)
{
    constexpr int num_grow = 3;

    std::swap(phi_old[lev], phi_new[lev]);

    MultiFab& S_new = phi_new[lev];

    const Real old_time = t_old[lev];
    const Real new_time = t_new[lev];
    const Real ctr_time = 0.5*(old_time+new_time);

    const auto dx = geom[lev].CellSizeArray();
    GpuArray<Real, AMREX_SPACEDIM> dtdx;
    for (int i=0; i<AMREX_SPACEDIM; ++i)
    {
        dtdx[i] = dt_lev/(dx[i]);
    }

    const Real* prob_lo = geom[lev].ProbLo();

    MultiFab fluxes[AMREX_SPACEDIM];
    if (do_reflux)
    {
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
        {
            BoxArray ba = grids[lev];
            ba.surroundingNodes(i);
            fluxes[i].define(ba, dmap[lev], S_new.nComp(), 0);
        }
    }

    // State with ghost cells
    MultiFab Sborder(grids[lev], dmap[lev], S_new.nComp(), num_grow);
    FillPatch(lev, time, Sborder, 0, Sborder.nComp());

    // Build temporary multiFabs to work on.
    Array<MultiFab, AMREX_SPACEDIM> fluxcalc;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        BoxArray ba = amrex::convert(S_new.boxArray(), IntVect::TheDimensionVector(idim));
        fluxcalc[idim].define (ba,S_new.DistributionMap(), S_new.nComp(), 0);
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
	for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{

        // ======== GET FACE VELOCITY =========
            GpuArray<Box, AMREX_SPACEDIM> nbx;
            AMREX_D_TERM(nbx[0] = mfi.nodaltilebox(0);,
                         nbx[1] = mfi.nodaltilebox(1);,
                         nbx[2] = mfi.nodaltilebox(2););

            AMREX_D_TERM(const Box& ngbxx = amrex::grow(mfi.nodaltilebox(0),1);,
                         const Box& ngbxy = amrex::grow(mfi.nodaltilebox(1),1);,
                         const Box& ngbxz = amrex::grow(mfi.nodaltilebox(2),1););

            GpuArray<Array4<Real>, AMREX_SPACEDIM> vel{ AMREX_D_DECL( facevel[lev][0].array(mfi),
                                                                      facevel[lev][1].array(mfi),
                                                                      facevel[lev][2].array(mfi)) };

        // ======== FLUX CALC AND UPDATE =========

	    const Box& bx = mfi.tilebox();
            const Box& gbx = amrex::grow(bx, 1);

            Array4<Real> statein  = Sborder.array(mfi);
            Array4<Real> stateout = S_new.array(mfi);

            GpuArray<Array4<Real>, AMREX_SPACEDIM> flux{ AMREX_D_DECL(fluxcalc[0].array(mfi),
                                                                      fluxcalc[1].array(mfi),
                                                                      fluxcalc[2].array(mfi)) };

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

            // compute new state (stateout) and scale fluxes based on face area.
            // ===========================

            // Do a conservative update 
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                conservative(i, j, k,
                             statein, stateout,
                             AMREX_D_DECL(flux[0], flux[1], flux[2]),
                             dtdx);
            });

            // Scale by face area in order to correctly reflux
            AMREX_D_TERM(
                         amrex::ParallelFor(amrex::growHi(bx, 0, 1),
                         [=] AMREX_GPU_DEVICE (int i, int j, int k)
                         {
                             flux_scale_x(i, j, k, flux[0], dt_lev, dx);
                         });,
 
                         amrex::ParallelFor(amrex::growHi(bx, 1, 1),
                         [=] AMREX_GPU_DEVICE (int i, int j, int k)
                         {
                             flux_scale_y(i, j, k, flux[1], dt_lev, dx);
                         });,

                         amrex::ParallelFor(amrex::growHi(bx, 2, 1),
                         [=] AMREX_GPU_DEVICE (int i, int j, int k)
                         {
                             flux_scale_z(i, j, k, flux[2], dt_lev, dx);
                         });
                        );

            if (do_reflux) {

                GpuArray<Array4<Real>, AMREX_SPACEDIM> fluxout{ AMREX_D_DECL(fluxes[0].array(mfi),
                                                                             fluxes[1].array(mfi),
                                                                             fluxes[2].array(mfi)) };

                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    amrex::ParallelFor(nbx[idim],
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        fluxout[idim](i,j,k) = flux[idim](i,j,k);
                    });
                }
            }
        }
    }

    // ======== CFL CHECK, MOVED OUTSIDE MFITER LOOP =========

    AMREX_D_TERM(Real umax = facevel[lev][0].norm0(0,0,false);,
                 Real vmax = facevel[lev][1].norm0(0,0,false);,
                 Real wmax = facevel[lev][2].norm0(0,0,false););

    if (AMREX_D_TERM(umax*dt_lev > dx[0], ||
                     vmax*dt_lev > dx[1], ||
                     wmax*dt_lev > dx[2]))
    {
#if (AMREX_SPACEDIM > 2)
        amrex::Print() << "umax = " << umax << ", vmax = " << vmax << ", wmax = " << wmax 
                       << ", dt = " << dt_lev << " dx = " << dx[0] << " " << dx[1] << " " << dx[2] << std::endl;
#else
        amrex::Print() << "umax = " << umax << ", vmax = " << vmax 
                       << ", dt = " << dt_lev << " dx = " << dx[0] << " " << dx[1] << " " << dx[2] << std::endl;
#endif
        amrex::Abort("CFL violation. use smaller adv.cfl.");
    }

    // ======== END OF GPU EDIT, (FOR NOW) =========

    // increment or decrement the flux registers by area and time-weighted fluxes
    // Note that the fluxes have already been scaled by dt and area
    // In this example we are solving phi_t = -div(+F)
    // The fluxes contain, e.g., F_{i+1/2,j} = (phi*u)_{i+1/2,j}
    // Keep this in mind when considering the different sign convention for updating
    // the flux registers from the coarse or fine grid perspective
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    if (do_reflux) { 
	if (flux_reg[lev+1]) {
	    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
	        // update the lev+1/lev flux register (index lev+1)   
	        flux_reg[lev+1]->CrseInit(fluxes[i],i,0,0,fluxes[i].nComp(), -1.0);
	    }	    
	}
	if (flux_reg[lev]) {
	    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
	        // update the lev/lev-1 flux register (index lev) 
		flux_reg[lev]->FineAdd(fluxes[i],i,0,0,fluxes[i].nComp(), 1.0);
	    }
	}
    }
}
