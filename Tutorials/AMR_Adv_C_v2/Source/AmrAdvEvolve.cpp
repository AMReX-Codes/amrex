
#include <limits>

#include <AmrAdv.H>
#include <AmrAdv_F.H>

void
AmrAdv::Evolve ()
{
    Real cur_time = t_new[0];

    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
	if (ParallelDescriptor::IOProcessor()) {
	    std::cout << "\nSTEP " << step+1 << " starts ..." << std::endl;
	}

	ComputeDt();

	int lev = 0;
	int iteration = 1;
	timeStep(lev, cur_time, iteration);

	cur_time += dt[0];

	if (ParallelDescriptor::IOProcessor()) {
	    std::cout << "STEP " << step+1 << " ends." << " TIME = " << cur_time << " DT = " << dt[0]
		      << std::endl;
	}

	// post coarsetimestep, io?

	if (cur_time >= stop_time - 1.e-6*dt[0]) break;
    }
}

void
AmrAdv::timeStep (int lev, Real time, int iteration)
{
    int lev_top = finest_level;

    for (int i = lev; i <= lev_top; ++i)
    {
	int old_finest = finest_level;

	// We may need to regrid
	{
	    
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
    isubstep[lev] = iteration;

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
    }

    // post time step stuff
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

    // State with ghost cells
    MultiFab Sborder(grids[lev], S_new.nComp(), num_grow, dmap[lev]);
    FillPatch(lev, time, Sborder, 0, Sborder.nComp());

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
			      D_DECL(BL_TO_FORTRAN(uface[0]),
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
