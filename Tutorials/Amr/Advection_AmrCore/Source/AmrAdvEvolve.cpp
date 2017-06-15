
#include <limits>

#include <AmrAdv.H>
#include <AmrAdv_F.H>

using namespace amrex;

void
AmrAdv::Evolve ()
{
    Real cur_time = t_new[0];
    int last_plot_file_step = 0;

    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
	if (ParallelDescriptor::IOProcessor()) {
	    std::cout << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;
	}

	ComputeDt();

	int lev = 0;
	int iteration = 1;
	timeStep(lev, cur_time, iteration);

	cur_time += dt[0];

	if (ParallelDescriptor::IOProcessor()) {
	    std::cout << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time << " DT = " << dt[0]
		      << std::endl;
	}

	// sync up time
	for (int lev = 0; lev <= finest_level; ++lev) {
	    t_new[lev] = cur_time;
	}

	if (plot_int > 0 && (step+1) % plot_int == 0) {
	    last_plot_file_step = step+1;
	    WritePlotFile();
	}

	if (cur_time >= stop_time - 1.e-6*dt[0]) break;
    }

    if (plot_int > 0 && istep[0] > last_plot_file_step) {
	WritePlotFile();
    }
}

void
AmrAdv::timeStep (int lev, Real time, int iteration)
{
    if (regrid_int > 0)  // We may need to regrid
    {
        static Array<int> last_regrid_step(max_level+1, 0);

        // regrid changes level "lev+1" so we don't regrid on max_level
        if (lev < max_level && istep[lev] > last_regrid_step[lev])
        {
            if (istep[lev] % regrid_int == 0)
            {
                // regrid could change finest_level if finest_level < max_level
                // so we save the previous finest level index
		int old_finest = finest_level; 
		regrid(lev, time);

		for (int k = lev; k <= finest_level; ++k) {
		    last_regrid_step[k] = istep[k];
		}

		for (int k = old_finest+1; k <= finest_level; ++k) {
		    dt[k] = dt[k-1] / MaxRefRatio(k-1);
		}
	    }
	}
    }

    if (Verbose() && ParallelDescriptor::IOProcessor()) {
	std::cout << "[Level " << lev 
		  << " step " << istep[lev]+1 << "] ";
	std::cout << "ADVANCE with dt = "
		  << dt[lev]
		  << std::endl;
    }

    Advance(lev, time, dt[lev], iteration, nsubsteps[lev]);

    ++istep[lev];

    if (Verbose() && ParallelDescriptor::IOProcessor())
    {
	std::cout << "[Level " << lev
		  << " step " << istep[lev] << "] ";
        std::cout << "Advanced "
                  << CountCells(lev)
                  << " cells"
                  << std::endl;
    }

    if (lev < finest_level)
    {
	for (int i = 1; i <= nsubsteps[lev+1]; ++i)
	{
	    timeStep(lev+1, time+(i-1)*dt[lev+1], i);
	}

	if (do_reflux)
	{
	    flux_reg[lev+1]->Reflux(*phi_new[lev], 1.0, 0, 0, phi_new[lev]->nComp(), geom[lev]);
	}

	AverageDownTo(lev); // average lev+1 down to lev
    }
    
}

void
AmrAdv::Advance (int lev, Real time, Real dt, int iteration, int ncycle)
{
    constexpr int num_grow = 3;

    std::swap(phi_old[lev], phi_new[lev]);
    t_old[lev] = t_new[lev];
    t_new[lev] += dt;

    MultiFab& S_new = *phi_new[lev];

    const Real old_time = t_old[lev];
    const Real new_time = t_new[lev];
    const Real ctr_time = 0.5*(old_time+new_time);

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

    MultiFab fluxes[BL_SPACEDIM];
    if (do_reflux)
    {
	for (int i = 0; i < BL_SPACEDIM; ++i)
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
#pragma omp parallel
#endif
    {
	FArrayBox flux[BL_SPACEDIM], uface[BL_SPACEDIM];

	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.tilebox();

	    const FArrayBox& statein = Sborder[mfi];
	    FArrayBox& stateout      =   S_new[mfi];

	    // Allocate fabs for fluxes and Godunov velocities.
	    for (int i = 0; i < BL_SPACEDIM ; i++) {
		const Box& bxtmp = amrex::surroundingNodes(bx,i);
		flux[i].resize(bxtmp,S_new.nComp());
		uface[i].resize(amrex::grow(bxtmp,1),1);
	    }

	    get_face_velocity(lev, ctr_time,
			      AMREX_D_DECL(BL_TO_FORTRAN(uface[0]),
				     BL_TO_FORTRAN(uface[1]),
				     BL_TO_FORTRAN(uface[2])),
			      dx, prob_lo);

            advect(time, bx.loVect(), bx.hiVect(),
		   BL_TO_FORTRAN_3D(statein), 
		   BL_TO_FORTRAN_3D(stateout),
		   AMREX_D_DECL(BL_TO_FORTRAN_3D(uface[0]),
			  BL_TO_FORTRAN_3D(uface[1]),
			  BL_TO_FORTRAN_3D(uface[2])),
		   AMREX_D_DECL(BL_TO_FORTRAN_3D(flux[0]), 
			  BL_TO_FORTRAN_3D(flux[1]), 
			  BL_TO_FORTRAN_3D(flux[2])), 
		   dx, dt);

	    if (do_reflux) {
		for (int i = 0; i < BL_SPACEDIM ; i++) {
		    fluxes[i][mfi].copy(flux[i],mfi.nodaltilebox(i));	  
		}
	    }
	}
    }

    if (do_reflux) { // Note that the fluxes have already been scaled by dt and area.
	if (flux_reg[lev]) {
	    for (int i = 0; i < BL_SPACEDIM; ++i) {
		flux_reg[lev]->FineAdd(fluxes[i],i,0,0,fluxes[i].nComp(), 1.0);
	    }
	}

	if (flux_reg[lev+1]) {
	    flux_reg[lev+1]->setVal(0.0);
	    for (int i = 0; i < BL_SPACEDIM; ++i) {
		flux_reg[lev+1]->CrseInit(fluxes[i],i,0,0,fluxes[i].nComp(), -1.0);
	    }	    
	}
    }
}

void
AmrAdv::ComputeDt ()
{
    Array<Real> dt_tmp(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
	dt_tmp[lev] = EstTimeStep(lev, true);
    }
    ParallelDescriptor::ReduceRealMin(&dt_tmp[0], dt_tmp.size());

    constexpr Real change_max = 1.1;
    Real dt_0 = dt_tmp[0];
    int n_factor = 1;
    for (int lev = 0; lev <= finest_level; ++lev) {
	dt_tmp[lev] = std::min(dt_tmp[lev], change_max*dt[lev]);
	n_factor *= nsubsteps[lev];
	dt_0 = std::min(dt_0, n_factor*dt_tmp[lev]);
    }

    // Limit dt's by the value of stop_time.
    const Real eps = 1.e-3*dt_0;
    if (t_new[0] + dt_0 > stop_time - eps) {
	dt_0 = stop_time - t_new[0];
    }

    dt[0] = dt_0;
    for (int lev = 1; lev <= finest_level; ++lev) {
	dt[lev] = dt[lev-1] / nsubsteps[lev];
    }
}

Real
AmrAdv::EstTimeStep (int lev, bool local) const
{
    BL_PROFILE("AmrAdv::EstTimeStep()");

    Real dt_est = std::numeric_limits<Real>::max();

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    const Real cur_time = t_new[lev];
    const MultiFab& S_new = *phi_new[lev];

#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_est)
#endif
    {
	FArrayBox uface[BL_SPACEDIM];

	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
	    for (int i = 0; i < BL_SPACEDIM ; i++) {
		const Box& bx = mfi.nodaltilebox(i);
		uface[i].resize(bx,1);
	    }

	    get_face_velocity(lev, cur_time,
			      AMREX_D_DECL(BL_TO_FORTRAN(uface[0]),
				     BL_TO_FORTRAN(uface[1]),
				     BL_TO_FORTRAN(uface[2])),
			      dx, prob_lo);

	    for (int i = 0; i < BL_SPACEDIM; ++i) {
		Real umax = uface[i].norm(0);
		if (umax > 1.e-100) {
		    dt_est = std::min(dt_est, dx[i] / umax);
		}
	    }
	}
    }

    if (!local) {
	ParallelDescriptor::ReduceRealMin(dt_est);
    }

    dt_est *= cfl;

    return dt_est;
}
