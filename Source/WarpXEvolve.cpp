
#include <cmath>
#include <limits>

#include <WarpX.H>
#include <WarpXConst.H>
#include <WarpX_f.H>

using namespace amrex;

void
WarpX::Evolve (int numsteps)
{
    BL_PROFILE("WarpX::Evolve()");

    Real cur_time = t_new[0];
    static int last_plot_file_step = 0;
    static int last_check_file_step = 0;

    int numsteps_max = (numsteps >= 0 && numsteps <= max_step) ? numsteps : max_step;
    bool max_time_reached = false;

    for (int step = istep[0]; step < numsteps_max && cur_time < stop_time; ++step)
    {
	if (ParallelDescriptor::IOProcessor()) {
	    std::cout << "\nSTEP " << step+1 << " starts ..." << std::endl;
	}

        if (ParallelDescriptor::NProcs() > 1)
           if (okToRegrid(step)) RegridBaseLevel();

	ComputeDt();

	// Advance level 0 by dt
	const int lev = 0;
	{
	    // At the beginning, we have B^{n-1/2} and E^{n}.
	    // Particles have p^{n-1/2} and x^{n}.

	    EvolveB(lev, 0.5*dt[lev]); // We now B^{n}

	    if (WarpX::nox > 1 || WarpX::noy > 1 || WarpX::noz > 1) {
		WarpX::FillBoundary(*Bfield[lev][0], geom[lev], Bx_nodal_flag);
		WarpX::FillBoundary(*Bfield[lev][1], geom[lev], By_nodal_flag);
		WarpX::FillBoundary(*Bfield[lev][2], geom[lev], Bz_nodal_flag);
		WarpX::FillBoundary(*Efield[lev][0], geom[lev], Ex_nodal_flag);
		WarpX::FillBoundary(*Efield[lev][1], geom[lev], Ey_nodal_flag);
		WarpX::FillBoundary(*Efield[lev][2], geom[lev], Ez_nodal_flag);
	    }

	    // Evolve particles to p^{n+1/2} and x^{n+1}
	    // Depose current, j^{n+1/2}
	    mypc->Evolve(lev,
			 *Efield[lev][0],*Efield[lev][1],*Efield[lev][2],
			 *Bfield[lev][0],*Bfield[lev][1],*Bfield[lev][2],
			 *current[lev][0],*current[lev][1],*current[lev][2], cur_time, dt[lev]);

	    mypc->Redistribute();  // Redistribute particles

	    EvolveB(lev, 0.5*dt[lev]); // We now B^{n+1/2}

	    // Fill B's ghost cells because of the next step of evolving E.
	    WarpX::FillBoundary(*Bfield[lev][0], geom[lev], Bx_nodal_flag);
	    WarpX::FillBoundary(*Bfield[lev][1], geom[lev], By_nodal_flag);
	    WarpX::FillBoundary(*Bfield[lev][2], geom[lev], Bz_nodal_flag);

	    EvolveE(lev, dt[lev]); // We now have E^{n+1}

	    ++istep[lev];
	}

	cur_time += dt[0];

	MoveWindow();

	if (ParallelDescriptor::IOProcessor()) {
	    std::cout << "STEP " << step+1 << " ends." << " TIME = " << cur_time << " DT = " << dt[0]
		      << std::endl;
	}

	// sync up time
	for (int i = 0; i <= finest_level; ++i) {
	    t_new[i] = cur_time;
	}

	if (plot_int > 0 && (step+1) % plot_int == 0) {
	    last_plot_file_step = step+1;
	    WritePlotFile();
	}

	if (check_int > 0 && (step+1) % check_int == 0) {
	    last_check_file_step = step+1;
	    WriteCheckPointFile();
	}

	if (cur_time >= stop_time - 1.e-6*dt[0]) {
	    max_time_reached = true;
	    break;
	}
    }

    if (plot_int > 0 && istep[0] > last_plot_file_step && (max_time_reached || istep[0] >= max_step)) {
	WritePlotFile();
    }

    if (check_int > 0 && istep[0] > last_check_file_step && (max_time_reached || istep[0] >= max_step)) {
	WriteCheckPointFile();
    }
}

void
WarpX::EvolveB (int lev, Real dt)
{
    BL_PROFILE("WarpX::EvolveB()");

    const Real* dx = geom[lev].CellSize();

    Real dtsdx[3];
#if (BL_SPACEDIM == 3)
    dtsdx[0] = dt / dx[0];
    dtsdx[1] = dt / dx[1];
    dtsdx[2] = dt / dx[2];
#elif (BL_SPACEDIM == 2)
    dtsdx[0] = dt / dx[0];
    dtsdx[1] = std::numeric_limits<Real>::quiet_NaN();
    dtsdx[2] = dt / dx[1];
#endif

    const int norder = 2;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(*Efield[lev][0],true); mfi.isValid(); ++mfi )
    {
	const Box& tbx = mfi.tilebox();
	WRPX_PXR_PUSH_BVEC(tbx.loVect(), tbx.hiVect(),
			   BL_TO_FORTRAN_3D((*Efield[lev][0])[mfi]),
			   BL_TO_FORTRAN_3D((*Efield[lev][1])[mfi]),
			   BL_TO_FORTRAN_3D((*Efield[lev][2])[mfi]),
			   BL_TO_FORTRAN_3D((*Bfield[lev][0])[mfi]),
			   BL_TO_FORTRAN_3D((*Bfield[lev][1])[mfi]),
			   BL_TO_FORTRAN_3D((*Bfield[lev][2])[mfi]),
			   dtsdx, dtsdx+1, dtsdx+2,
			   &norder);
    }
}

void
WarpX::EvolveE (int lev, Real dt)
{
    BL_PROFILE("WarpX::EvolveE()");

    Real mu_c2_dt = (PhysConst::mu0*PhysConst::c*PhysConst::c) * dt;

    const Real* dx = geom[lev].CellSize();

    Real dtsdx_c2[3];
#if (BL_SPACEDIM == 3)
    dtsdx_c2[0] = (PhysConst::c*PhysConst::c) * dt / dx[0];
    dtsdx_c2[1] = (PhysConst::c*PhysConst::c) * dt / dx[1];
    dtsdx_c2[2] = (PhysConst::c*PhysConst::c) * dt / dx[2];
#else
    dtsdx_c2[0] = (PhysConst::c*PhysConst::c) * dt / dx[0];
    dtsdx_c2[1] = std::numeric_limits<Real>::quiet_NaN();
    dtsdx_c2[2] = (PhysConst::c*PhysConst::c) * dt / dx[1];
#endif

    const int norder = 2;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(*Efield[lev][0],true); mfi.isValid(); ++mfi )
    {
	const Box& tbx = mfi.tilebox();
	WRPX_PXR_PUSH_EVEC(tbx.loVect(), tbx.hiVect(),
			   BL_TO_FORTRAN_3D((*Efield[lev][0])[mfi]),
			   BL_TO_FORTRAN_3D((*Efield[lev][1])[mfi]),
			   BL_TO_FORTRAN_3D((*Efield[lev][2])[mfi]),
			   BL_TO_FORTRAN_3D((*Bfield[lev][0])[mfi]),
			   BL_TO_FORTRAN_3D((*Bfield[lev][1])[mfi]),
			   BL_TO_FORTRAN_3D((*Bfield[lev][2])[mfi]),
			   BL_TO_FORTRAN_3D((*current[lev][0])[mfi]),
			   BL_TO_FORTRAN_3D((*current[lev][1])[mfi]),
			   BL_TO_FORTRAN_3D((*current[lev][2])[mfi]),
			   &mu_c2_dt,
			   dtsdx_c2, dtsdx_c2+1, dtsdx_c2+2,
			   &norder);
    }
}

void
WarpX::ComputeDt ()
{
    Array<Real> dt_tmp(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
	const Real* dx = geom[lev].CellSize();
	dt_tmp[lev]  = cfl * 1./( std::sqrt(D_TERM(  1./(dx[0]*dx[0]),
						   + 1./(dx[1]*dx[1]),
						   + 1./(dx[2]*dx[2]))) * PhysConst::c );
    }

    // Limit dt's by the value of stop_time.
    Real dt_0 = dt_tmp[0];
    const Real eps = 1.e-3*dt_0;
    if (t_new[0] + dt_0 > stop_time - eps) {
	dt_0 = stop_time - t_new[0];
    }

    dt[0] = dt_0;
    for (int lev = 1; lev <= finest_level; ++lev) {
	dt[lev] = dt[lev-1] / nsubsteps[lev];
    }
}

void
WarpX::InjectPlasma (int num_shift, int dir)
{
    if(do_plasma_injection)
    {
        // particleBox encloses the cells where we generate particles
        Box particleBox = geom[0].Domain();
        int domainLength = particleBox.length(dir);
        int sign = (num_shift < 0) ? -1 : 1;
        particleBox.shift(dir, sign*(domainLength - std::abs(num_shift)));
        particleBox &= geom[0].Domain();

        const Real* dx  = geom[0].CellSize();

        const int lev = 0;

        std::array<Real,PIdx::nattribs> attribs;
        attribs.fill(0.0);

        WarpXParticleContainer& pc0 = mypc->GetParticleContainer(0);        
        for (MFIter mfi = pc0.MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            const Box& tile_box     = mfi.tilebox();
            const Box& intersectBox = tile_box & particleBox;
            if (intersectBox.ok())
            {
                RealBox real_box { intersectBox, dx, geom[0].ProbLo() };

                const int grid_id = mfi.index();
                const int tile_id = mfi.LocalTileIndex();

                for (int ispec=0; ispec < num_injected_species; ispec++)
                {
                    int ispecies = injected_plasma_species[ispec];
                    WarpXParticleContainer& pc = mypc->GetParticleContainer(ispecies);

                    const IntVect& boxlo = tile_box.smallEnd();
                    for (IntVect iv = boxlo; iv <= tile_box.bigEnd(); tile_box.next(iv))
                    {
                        for (int i_part=0; i_part < injected_plasma_ppc[ispec]; i_part++)
                        {
                            Real particle_shift = (0.5+i_part)/injected_plasma_ppc[ispec];
#if (BL_SPACEDIM == 3)
                            Real x = real_box.lo(0) + (iv[0]-boxlo[0] + particle_shift)*dx[0];
                            Real y = real_box.lo(1) + (iv[1]-boxlo[1] + particle_shift)*dx[1];
                            Real z = real_box.lo(2) + (iv[2]-boxlo[2] + particle_shift)*dx[2];
#elif (BL_SPACEDIM == 2)
                            Real x = real_box.lo(0) + (iv[0]-boxlo[0] + particle_shift)*dx[0];
                            Real y = 0.0;
                            Real z = real_box.lo(1) + (iv[1]-boxlo[1] + particle_shift)*dx[1];
#endif

                            Real weight = injected_plasma_density[ispec];
#if BL_SPACEDIM==3
                            weight *= dx[0]*dx[1]*dx[2]/injected_plasma_ppc[ispec];
#elif BL_SPACEDIM==2
                            weight *= dx[0]*dx[1]/injected_plasma_ppc[ispec];
#endif
                            attribs[PIdx::w] = weight;

                            pc.AddOneParticle(lev, grid_id, tile_id, x, y, z, attribs);
                        }
                    }
                }
            }
        }
    }
}

void
WarpX::MoveWindow ()
{

  if (do_moving_window == 0) return;

  // compute the number of cells to shift
  int dir = moving_window_dir;
  Real new_lo[BL_SPACEDIM];
  Real new_hi[BL_SPACEDIM];
  const Real* current_lo = geom[0].ProbLo();
  const Real* current_hi = geom[0].ProbHi();
  const Real* dx = geom[0].CellSize();
  moving_window_x += moving_window_v * dt[0];
  int num_shift = (moving_window_x - current_lo[dir]) / dx[dir];

  if (num_shift == 0) return;

  // update the problem domain
  for (int i=0; i<BL_SPACEDIM; i++) {
    new_lo[i] = current_lo[i];
    new_hi[i] = current_hi[i];
  }
  new_lo[dir] = current_lo[dir] + num_shift * dx[dir];
  new_hi[dir] = current_hi[dir] + num_shift * dx[dir];
  RealBox new_box(new_lo, new_hi);
  geom[0].ProbDomain(new_box);

  // shift the mesh fields (Note - only on level 0 for now)
  shiftMF(*Bfield[0][0], geom[0], num_shift, dir, Bx_nodal_flag);
  shiftMF(*Bfield[0][1], geom[0], num_shift, dir, By_nodal_flag);
  shiftMF(*Bfield[0][2], geom[0], num_shift, dir, Bz_nodal_flag);
  shiftMF(*Efield[0][0], geom[0], num_shift, dir, Ex_nodal_flag);
  shiftMF(*Efield[0][1], geom[0], num_shift, dir, Ey_nodal_flag);
  shiftMF(*Efield[0][2], geom[0], num_shift, dir, Ez_nodal_flag);

  InjectPlasma(num_shift, dir);

  // Redistribute (note - this removes particles that are outside of the box)
  mypc->Redistribute();
}
