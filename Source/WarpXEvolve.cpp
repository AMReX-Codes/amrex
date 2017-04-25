
#include <cmath>
#include <limits>

#include <WarpX.H>
#include <WarpXConst.H>
#include <WarpX_f.H>
#include <WarpX_py.H>

using namespace amrex;

void
WarpX::Evolve (int numsteps)
{
    BL_PROFILE("WarpX::Evolve()");

    Real cur_time = t_new[0];
    static int last_plot_file_step = 0;
    static int last_check_file_step = 0;

    int numsteps_max;
    if (numsteps < 0) {  // Note that the default argument is numsteps = -1
        numsteps_max = max_step;
    } else {
        numsteps_max = std::min(istep[0]+numsteps, max_step);
    }

    bool max_time_reached = false;
    bool last_step = false;

    for (int step = istep[0]; step < numsteps_max && cur_time < stop_time; ++step)
    {
        if (warpx_py_print_step) {
            warpx_py_print_step(step);
        }

	// Start loop on time steps
        amrex::Print() << "\nSTEP " << step+1 << " starts ...\n";

        if (ParallelDescriptor::NProcs() > 1) {
           if (okToRegrid(step)) RegridBaseLevel();
        }

	ComputeDt();

	// Advance level 0 by dt
	const int lev = 0;
	{
	    // At the beginning, we have B^{n-1/2} and E^{n-1/2}.
	    // Particles have p^{n-1/2} and x^{n-1/2}.

	    // Beyond one step, we have B^{n-1/2} and E^{n}.
	    // Particles have p^{n-1/2} and x^{n}.

	    if (is_synchronized) {
	        // on first step, push E and X by 0.5*dt
                WarpX::FillBoundaryB( lev, true );
	        EvolveE(lev, 0.5*dt[lev]);
	        mypc->PushX(lev, 0.5*dt[lev]);
                mypc->Redistribute();  // Redistribute particles
                is_synchronized = false;
	    }

	    EvolveB(lev, 0.5*dt[lev]); // We now B^{n}

            WarpX::FillBoundaryB( lev, false );
            WarpX::FillBoundaryE( lev, false );

	    // Evolve particles to p^{n+1/2} and x^{n+1}
	    // Depose current, j^{n+1/2}
	    mypc->Evolve(lev,
			 *Efield[lev][0],*Efield[lev][1],*Efield[lev][2],
			 *Bfield[lev][0],*Bfield[lev][1],*Bfield[lev][2],
			 *current[lev][0],*current[lev][1],*current[lev][2], cur_time, dt[lev]);

	    EvolveB(lev, 0.5*dt[lev]); // We now B^{n+1/2}

	    // Fill B's ghost cells because of the next step of evolving E.
            WarpX::FillBoundaryB( lev, true );

   	    if (cur_time + dt[0] >= stop_time - 1.e-6*dt[0] || step == numsteps_max-1) {
   	        // on last step, push by only 0.5*dt to synchronize all at n+1/2
	        EvolveE(lev, 0.5*dt[lev]); // We now have E^{n+1/2}
	        mypc->PushX(lev, -0.5*dt[lev]);
                is_synchronized = true;
            } else {
	        EvolveE(lev, dt[lev]); // We now have E^{n+1}
	    }

	    mypc->Redistribute();  // Redistribute particles

	    ++istep[lev];
	}

	cur_time += dt[0];

        bool to_make_plot = (plot_int > 0) && ((step+1) % plot_int == 0);

        bool move_j = is_synchronized || to_make_plot;
        // If is_synchronized we need to shift j too so that next step we can evolve E by dt/2.
        // We might need to move j because we are going to make a plotfile.
	MoveWindow(move_j);

        amrex::Print()<< "STEP " << step+1 << " ends." << " TIME = " << cur_time
                      << " DT = " << dt[0] << "\n";

	// sync up time
	for (int i = 0; i <= finest_level; ++i) {
	    t_new[i] = cur_time;
	}

	if (to_make_plot) {
            WarpX::FillBoundaryB( lev, false );
            WarpX::FillBoundaryE( lev, false );
            mypc->FieldGather(lev,
                              *Efield[lev][0],*Efield[lev][1],*Efield[lev][2],
                              *Bfield[lev][0],*Bfield[lev][1],*Bfield[lev][2]);
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

	// End loop on time steps
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

    // Parameters of the solver: order and mesh spacing
    const int norder = 2;

    const std::array<Real,3>& dx = WarpX::CellSize(lev);
    const std::array<Real,3> dtsdx {dt/dx[0], dt/dx[1], dt/dx[2]};

    // Loop through the grids, and over the tiles within each grid
#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(*Bfield[lev][0],true); mfi.isValid(); ++mfi )
    {
        const Box& tbx  = mfi.tilebox(Bx_nodal_flag);
        const Box& tby  = mfi.tilebox(By_nodal_flag);
        const Box& tbz  = mfi.tilebox(Bz_nodal_flag);

	// Call picsar routine for each tile
	WRPX_PXR_PUSH_BVEC(
	    tbx.loVect(), tbx.hiVect(),
	    tby.loVect(), tby.hiVect(),
	    tbz.loVect(), tbz.hiVect(),
	    BL_TO_FORTRAN_3D((*Efield[lev][0])[mfi]),
	    BL_TO_FORTRAN_3D((*Efield[lev][1])[mfi]),
	    BL_TO_FORTRAN_3D((*Efield[lev][2])[mfi]),
	    BL_TO_FORTRAN_3D((*Bfield[lev][0])[mfi]),
	    BL_TO_FORTRAN_3D((*Bfield[lev][1])[mfi]),
	    BL_TO_FORTRAN_3D((*Bfield[lev][2])[mfi]),
	    &dtsdx[0], &dtsdx[1], &dtsdx[2],
	    &norder);
    }

    if (do_pml && lev == 0)
    {
#if 0
        const Geometry& gm = Geom(0);
        const Box& domain = gm.Domain();
        const int* dlo = domain.loVect();
        const int* dhi = domain.hiVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(*pmlB[0]); mfi.isValid(); ++mfi )
        {
            const Box& ccbx = amrex::enclosedCells(mfi.validbox());
            const int* cclo = ccbx.loVect();
            const int* cchi = ccbx.hiVect();
            Box Bx_bx = amrex::convert(ccbx, Bx_nodal_flag); 
            Box By_bx = amrex::convert(ccbx, By_nodal_flag); 
            Box Bz_bx = amrex::convert(ccbx, Bz_nodal_flag);
        }
#endif
    }
}

void
WarpX::EvolveE (int lev, Real dt)
{
    BL_PROFILE("WarpX::EvolveE()");

    // Parameters of the solver: order and mesh spacing
    const int norder = 2;

    const std::array<Real,3>& dx = WarpX::CellSize(lev);
    const Real mu_c2_dt = (PhysConst::mu0*PhysConst::c*PhysConst::c) * dt;
    const Real foo = (PhysConst::c*PhysConst::c) * dt;
    const std::array<Real,3> dtsdx_c2 {foo/dx[0], foo/dx[1], foo/dx[2]};

  // Loop through the grids, and over the tiles within each grid
#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(*Efield[lev][0],true); mfi.isValid(); ++mfi )
    {
	const Box& tex  = mfi.tilebox(Ex_nodal_flag);
	const Box& tey  = mfi.tilebox(Ey_nodal_flag);
	const Box& tez  = mfi.tilebox(Ez_nodal_flag);

  // Call picsar routine for each tile
	WRPX_PXR_PUSH_EVEC(
	    tex.loVect(), tex.hiVect(),
	    tey.loVect(), tey.hiVect(),
	    tez.loVect(), tez.hiVect(),
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
	    &dtsdx_c2[0], &dtsdx_c2[1], &dtsdx_c2[2],
	    &norder);
    }
}

void
WarpX::FillBoundaryE(int lev, bool force)
{
    if (force || WarpX::nox > 1 || WarpX::noy > 1 || WarpX::noz > 1)
    {
        if (do_pml && lev == 0) {
            WarpX::ExchangeWithPML(*Efield[lev][0], *pml_E[0], geom[lev]);
            WarpX::ExchangeWithPML(*Efield[lev][1], *pml_E[1], geom[lev]);
            WarpX::ExchangeWithPML(*Efield[lev][2], *pml_E[2], geom[lev]);

            (*pml_E[0]).FillBoundary( geom[lev].periodicity() );
            (*pml_E[1]).FillBoundary( geom[lev].periodicity() );
            (*pml_E[2]).FillBoundary( geom[lev].periodicity() );
        }

        (*Efield[lev][0]).FillBoundary( geom[lev].periodicity() );
        (*Efield[lev][1]).FillBoundary( geom[lev].periodicity() );
        (*Efield[lev][2]).FillBoundary( geom[lev].periodicity() );
    }
}

void
WarpX::FillBoundaryB(int lev, bool force)
{
    if (force || WarpX::nox > 1 || WarpX::noy > 1 || WarpX::noz > 1)
    {
        if (do_pml && lev == 0) {
            WarpX::ExchangeWithPML(*Bfield[lev][0], *pml_B[0], geom[lev]);
            WarpX::ExchangeWithPML(*Bfield[lev][1], *pml_B[1], geom[lev]);
            WarpX::ExchangeWithPML(*Bfield[lev][2], *pml_B[2], geom[lev]);

            (*pml_B[0]).FillBoundary( geom[lev].periodicity() );
            (*pml_B[1]).FillBoundary( geom[lev].periodicity() );
            (*pml_B[2]).FillBoundary( geom[lev].periodicity() );
        }

        (*Bfield[lev][0]).FillBoundary( geom[lev].periodicity() );
        (*Bfield[lev][1]).FillBoundary( geom[lev].periodicity() );
        (*Bfield[lev][2]).FillBoundary( geom[lev].periodicity() );
    }
}

void
WarpX::ExchangeWithPML (MultiFab& regmf, MultiFab& pmlmf, const Geometry& gm)
{
    if (!do_pml) return;

    // Copy from PML to regular data

    MultiFab totpmlmf(pmlmf.boxArray(), pmlmf.DistributionMap(), 1, 0);
    MultiFab::LinComb(totpmlmf, 1.0, pmlmf, 0, 1.0, pmlmf, 1, 0, 1, 0);

    MultiFab tmpregmf(regmf.boxArray(), regmf.DistributionMap(), 1, regmf.nGrow());
    tmpregmf.copy(totpmlmf, 0, 0, 1, 0, regmf.nGrow(), gm.periodicity());

    Box dom = gm.Domain();
    dom.convert(regmf.ixType());
    for (int idim=0; idim < BL_SPACEDIM; ++idim) {
        if (Geometry::isPeriodic(idim)) {
            dom.grow(idim, regmf.nGrow());
        }
    }
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(regmf); mfi.isValid(); ++mfi)
    {
        const BoxList& bl = amrex::boxDiff(mfi.fabbox(), dom);
        for (const Box& bx : bl)
        {
            regmf[mfi].copy(tmpregmf[mfi], bx, 0, bx, 0, 1);
        }
    }

    // Copy from regular data to PML's first component
    pmlmf.copy (regmf, 0, 0, 1, 0, pmlmf.nGrow(), gm.periodicity());
    // Zero out the second component
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(pmlmf); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.fabbox();
        bx &= dom;
        if (bx.ok()) {
            pmlmf[mfi].setVal(0.0, bx, 1, 1);
        }
    }
}

void
WarpX::PushParticlesandDepose(int lev, Real cur_time)
{
    // Evolve particles to p^{n+1/2} and x^{n+1}
    // Depose current, j^{n+1/2}
    mypc->Evolve(lev,
		 *Efield[lev][0],*Efield[lev][1],*Efield[lev][2],
		 *Bfield[lev][0],*Bfield[lev][1],*Bfield[lev][2],
		 *current[lev][0],*current[lev][1],*current[lev][2], cur_time, dt[lev]);
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
        const int lev = 0;

        // particleBox encloses the cells where we generate particles
        Box particleBox = geom[lev].Domain();
        int domainLength = particleBox.length(dir);
        int sign = (num_shift < 0) ? -1 : 1;
        particleBox.shift(dir, sign*(domainLength - std::abs(num_shift)));
        particleBox &= geom[lev].Domain();

        const Real* dx = geom[lev].CellSize();

        for (int i = 0; i < num_injected_species; ++i) {
            int ispecies = injected_plasma_species[i];
            WarpXParticleContainer& pc = mypc->GetParticleContainer(ispecies);
            auto& ppc = dynamic_cast<PhysicalParticleContainer&>(pc);
            ppc.AddParticles(lev, particleBox);
        }
    }
}

void
WarpX::MoveWindow (bool move_j)
{

    if (do_moving_window == 0) return;

    const int lev = 0;

    // compute the number of cells to shift
    int dir = moving_window_dir;
    Real new_lo[BL_SPACEDIM];
    Real new_hi[BL_SPACEDIM];
    const Real* current_lo = geom[lev].ProbLo();
    const Real* current_hi = geom[lev].ProbHi();
    const Real* dx = geom[lev].CellSize();
    moving_window_x += moving_window_v * dt[lev];
    int num_shift = static_cast<int>((moving_window_x - current_lo[dir]) / dx[dir]);

    if (num_shift == 0) return;

    // update the problem domain
    for (int i=0; i<BL_SPACEDIM; i++) {
        new_lo[i] = current_lo[i];
        new_hi[i] = current_hi[i];
    }
    new_lo[dir] = current_lo[dir] + num_shift * dx[dir];
    new_hi[dir] = current_hi[dir] + num_shift * dx[dir];
    RealBox new_box(new_lo, new_hi);
    geom[lev].ProbDomain(new_box);

    // shift the mesh fields (Note - only on level 0 for now)
    shiftMF(*Bfield[lev][0], geom[lev], num_shift, dir);
    shiftMF(*Bfield[lev][1], geom[lev], num_shift, dir);
    shiftMF(*Bfield[lev][2], geom[lev], num_shift, dir);
    shiftMF(*Efield[lev][0], geom[lev], num_shift, dir);
    shiftMF(*Efield[lev][1], geom[lev], num_shift, dir);
    shiftMF(*Efield[lev][2], geom[lev], num_shift, dir);
    if (move_j) {
        shiftMF(*current[lev][0], geom[lev], num_shift, dir);
        shiftMF(*current[lev][1], geom[lev], num_shift, dir);
        shiftMF(*current[lev][2], geom[lev], num_shift, dir);
    }

    InjectPlasma(num_shift, dir);

    // Redistribute (note - this removes particles that are outside of the box)
    mypc->Redistribute();
}
