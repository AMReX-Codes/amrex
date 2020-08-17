#include <AmrCoreAdv.H>
#include <Kernels.H>

using namespace amrex;

// Advance a single level for a single time step, updates flux registers
void
AmrCoreAdv::AdvancePhiAtLevel (int lev, Real time, Real dt_lev, int /*iteration*/, int /*ncycle*/)
{
    constexpr int num_grow = 3;

    std::swap(phi_old[lev], phi_new[lev]);

    MultiFab& S_new = phi_new[lev];

    const Real dx = geom[lev].CellSize(0);
    const Real dy = geom[lev].CellSize(1);
    const Real dz = (AMREX_SPACEDIM == 2) ? Real(1.0) : geom[lev].CellSize(2);
    AMREX_D_TERM(Real dtdx = dt_lev/dx;,
                 Real dtdy = dt_lev/dy;,
                 Real dtdz = dt_lev/dz);

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

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FArrayBox tmpfab;
	for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
            AMREX_ASSERT(S_new.nComp() == 1);

        // ======== GET FACE VELOCITY =========

            AMREX_D_TERM(Array4<Real const> velx = facevel[lev][0].const_array(mfi);,
                         Array4<Real const> vely = facevel[lev][1].const_array(mfi);,
                         Array4<Real const> velz = facevel[lev][2].const_array(mfi));

        // ======== FLUX CALC AND UPDATE =========

	    const Box& bx = mfi.tilebox();
            const Box& gbx = amrex::grow(bx, 1);

            Array4<Real const> statein  = Sborder.const_array(mfi);
            Array4<Real      > stateout = S_new.array(mfi);

            int ntmpcomps = (AMREX_SPACEDIM == 2) ? 6 : 14;
            tmpfab.resize(amrex::grow(bx,2),ntmpcomps);
            Elixir tmpeli = tmpfab.elixir();
            int itmp = 0;

            AMREX_D_TERM(Array4<Real> tfluxx = tmpfab.array(itmp++);,
                         Array4<Real> tfluxy = tmpfab.array(itmp++);,
                         Array4<Real> tfluxz = tmpfab.array(itmp++));

            Array4<Real> slope2 = tmpfab.array(itmp++);
            Array4<Real const> slope2_c = slope2;
            Array4<Real> slope4 = tmpfab.array(itmp++);
            Array4<Real const> slope4_c = slope4;

            // compute longitudinal fluxes
            // ===========================

            // x -------------------------
            Array4<Real> phix = tmpfab.array(itmp++);
            Array4<Real const> phix_c = phix;

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
            Array4<Real> phiy = tmpfab.array(itmp++);
            Array4<Real const> phiy_c = phiy;

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
            Array4<Real> phiz = tmpfab.array(itmp++);
            Array4<Real const> phiz_c = phiz;

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
            Array4<Real> phix_y = tmpfab.array(itmp++);
            Array4<Real const> phix_y_c = phix_y;

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
            Array4<Real> phix_z = tmpfab.array(itmp++);
            Array4<Real const> phix_z_c = phix_z;

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
            Array4<Real> phiy_x = tmpfab.array(itmp++);
            Array4<Real const> phiy_x_c = phiy_x;

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
            Array4<Real> phiy_z = tmpfab.array(itmp++);
            Array4<Real const> phiy_z_c = phiy_z;

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
            Array4<Real> phiz_x = tmpfab.array(itmp++);
            Array4<Real const> phiz_x_c = phiz_x;

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
            Array4<Real> phiz_y = tmpfab.array(itmp++);
            Array4<Real const> phiz_y_c = phiz_y;

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
            amrex::ParallelFor(amrex::surroundingNodes(bx,Direction::x),
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                create_flux_x(i, j, k, tfluxx,
                              AMREX_D_DECL(velx,vely,velz),
#if (AMREX_SPACEDIM == 3)
                              phix_c, phiy_z_c, phiz_y_c,
                              dtdy, dtdz);
#else
                              phix_c, phiy_c,
                              dtdy);
#endif
            });

            amrex::ParallelFor(amrex::surroundingNodes(bx,Direction::y),
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                create_flux_y(i, j, k, tfluxy,
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
            amrex::ParallelFor(amrex::surroundingNodes(bx,Direction::z),
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                create_flux_z(i, j, k, tfluxz,
                              velx, vely, velz,
                              phiz_c, phix_y_c, phiy_x_c,
                              dtdx, dtdy);
            });
#endif
            AMREX_ASSERT(itmp == ntmpcomps);

            // compute new state (stateout) and scale fluxes based on face area.
            // ===========================

            AMREX_D_TERM(Array4<Real const> tfluxx_c = tfluxx;,
                         Array4<Real const> tfluxy_c = tfluxy;,
                         Array4<Real const> tfluxz_c = tfluxz);

            // Do a conservative update 
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                conservative(i, j, k,
                             stateout, statein,
                             AMREX_D_DECL(tfluxx_c,tfluxy_c,tfluxz_c),
                             AMREX_D_DECL(dtdx,dtdy,dtdz));
            });

            if (do_reflux)
            {
                // Scale by face area in order to correctly reflux
                amrex::ParallelFor(
                    AMREX_D_DECL(amrex::surroundingNodes(bx,Direction::x),
                                 amrex::surroundingNodes(bx,Direction::y),
                                 amrex::surroundingNodes(bx,Direction::z)),
                    AMREX_D_DECL([=] AMREX_GPU_DEVICE (int i, int j, int k)
                                 {
                                     tfluxx(i,j,k) *= dt_lev*dy*dz;
                                 },
                                 [=] AMREX_GPU_DEVICE (int i, int j, int k)
                                 {
                                     tfluxy(i,j,k) *= dt_lev*dx*dz;
                                 },
                                 [=] AMREX_GPU_DEVICE (int i, int j, int k)
                                 {
                                     tfluxz(i,j,k) *= dt_lev*dx*dy;
                                 }));
 
                // Copy into Flux MultiFab
                AMREX_D_TERM(Array4<Real> fluxx = fluxes[0].array(mfi);,
                             Array4<Real> fluxy = fluxes[1].array(mfi);,
                             Array4<Real> fluxz = fluxes[2].array(mfi));
                amrex::ParallelFor(
                    AMREX_D_DECL(mfi.nodaltilebox(0),
                                 mfi.nodaltilebox(1),
                                 mfi.nodaltilebox(2)),
                    AMREX_D_DECL([=] AMREX_GPU_DEVICE (int i, int j, int k)
                                 {
                                     fluxx(i,j,k) = tfluxx_c(i,j,k);
                                 },
                                 [=] AMREX_GPU_DEVICE (int i, int j, int k)
                                 {
                                     fluxy(i,j,k) = tfluxy_c(i,j,k);
                                 },
                                 [=] AMREX_GPU_DEVICE (int i, int j, int k)
                                 {
                                     fluxz(i,j,k) = tfluxz_c(i,j,k);
                                 }));
            }
        }
    }

    // ======== CFL CHECK, MOVED OUTSIDE MFITER LOOP =========

    AMREX_D_TERM(Real umax = facevel[lev][0].norminf(0,0,true);,
                 Real vmax = facevel[lev][1].norminf(0,0,true);,
                 Real wmax = facevel[lev][2].norminf(0,0,true));

    if (AMREX_D_TERM(umax*dt_lev > dx, ||
                     vmax*dt_lev > dy, ||
                     wmax*dt_lev > dz))
    {
#if (AMREX_SPACEDIM > 2)
        amrex::AllPrint() << "umax = " << umax << ", vmax = " << vmax << ", wmax = " << wmax
                          << ", dt = " << dt_lev << " dx = " << dx << " " << dy << " " << dz << std::endl;
#else
        amrex::AllPrint() << "umax = " << umax << ", vmax = " << vmax
                          << ", dt = " << dt_lev << " dx = " << dx << " " << dy << " " << dz << std::endl;
#endif
        amrex::Abort("CFL violation. use smaller adv.cfl.");
    }

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
