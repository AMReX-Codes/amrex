
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>

#include <AmrAdv.H>
#include <AmrAdvBC.H>
#include <AmrAdv_F.H>

using namespace amrex;

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
AmrAdv::AmrAdv ()
{
    ReadParameters();

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = maxLevel() + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= maxLevel(); ++lev) {
	nsubsteps[lev] = MaxRefRatio(lev-1);
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    phi_new.resize(nlevs_max);
    phi_old.resize(nlevs_max);

    flux_reg.resize(nlevs_max+1);
}

AmrAdv::~AmrAdv ()
{

}

// advance solution to final time
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

// initializes multilevel data
void
AmrAdv::InitData ()
{
    const Real time = 0.0;
    InitFromScratch(time);
    AverageDown();

    if (plot_int > 0) {
        WritePlotFile();
    }
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void
AmrAdv::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
				const DistributionMapping& dm)
{
    const int ncomp = phi_new[lev-1]->nComp();
    const int nghost = phi_new[lev-1]->nGrow();
    
    phi_new[lev].reset(new MultiFab(ba, dm, ncomp, nghost));
    phi_old[lev].reset(new MultiFab(ba, dm, ncomp, nghost));

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }

    FillCoarsePatch(lev, time, *phi_new[lev], 0, ncomp);
}

// Remake an existing level using provided BoxArray and DistributionMapping and 
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
AmrAdv::RemakeLevel (int lev, Real time, const BoxArray& ba,
		     const DistributionMapping& dm)
{
    const int ncomp = phi_new[lev]->nComp();
    const int nghost = phi_new[lev]->nGrow();

#if __cplusplus >= 201402L
    auto new_state = std::make_unique<MultiFab>(ba, dm, ncomp, nghost);
    auto old_state = std::make_unique<MultiFab>(ba, dm, ncomp, nghost);
#else
    std::unique_ptr<MultiFab> new_state(new MultiFab(ba, dm, ncomp, nghost));
    std::unique_ptr<MultiFab> old_state(new MultiFab(ba, dm, ncomp, nghost));
#endif

    FillPatch(lev, time, *new_state, 0, ncomp);

    std::swap(new_state, phi_new[lev]);
    std::swap(old_state, phi_old[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }    
}

// Delete level data
// overrides the pure virtual function in AmrCore
void
AmrAdv::ClearLevel (int lev)
{
    phi_new[lev].reset(nullptr);
    phi_old[lev].reset(nullptr);
    flux_reg[lev].reset(nullptr);
}

// initialize data using a fortran routine to compute initial state
// overrides the pure virtual function in AmrCore
void AmrAdv::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
				      const DistributionMapping& dm)
{
    const int ncomp = 1;
    const int nghost = 0;

    phi_new[lev].reset(new MultiFab(ba, dm, ncomp, nghost));
    phi_old[lev].reset(new MultiFab(ba, dm, ncomp, nghost));

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    Real cur_time = t_new[lev];

    MultiFab& state = *phi_new[lev];

    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const int* lo  = box.loVect();
        const int* hi  = box.hiVect();

	initdata(lev, cur_time, ARLIM_3D(lo), ARLIM_3D(hi),
		 BL_TO_FORTRAN_3D(state[mfi]), ZFILL(dx),
		 ZFILL(prob_lo));
    }
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
AmrAdv::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
    static bool first = true;
    static Array<Real> phierr;

    if (first)
    {
	first = false;
	ParmParse pp("adv");
	int n = pp.countval("phierr");
	if (n > 0) {
	    pp.getarr("phierr", phierr, 0, n);
	}
    }

    if (lev >= phierr.size()) return;

    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

    const MultiFab& state = *phi_new[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<int>  itags;
	
	for (MFIter mfi(state,true); mfi.isValid(); ++mfi)
	{
	    const Box&  tilebx  = mfi.tilebox();

            TagBox&     tagfab  = tags[mfi];
	    
	    // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
	    // So we are going to get a temporary integer array.
	    tagfab.get_itags(itags, tilebx);
	    
            // data pointer and index space
	    int*        tptr    = itags.dataPtr();
	    const int*  tlo     = tilebx.loVect();
	    const int*  thi     = tilebx.hiVect();

	    state_error(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
			BL_TO_FORTRAN_3D(state[mfi]),
			&tagval, &clearval, 
			ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()), 
			ZFILL(dx), ZFILL(prob_lo), &time, &phierr[lev]);
	    //
	    // Now update the tags in the TagBox.
	    //
	    tagfab.tags_and_untags(itags, tilebx);
	}
    }
}

// read in some parameters from inputs file
void
AmrAdv::ReadParameters ()
{
    {
	ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
	pp.query("max_step", max_step);
	pp.query("stop_time", stop_time);
    }

    {
	ParmParse pp("amr"); // Traditionally, these have prefix, amr.

	pp.query("regrid_int", regrid_int);
	pp.query("plot_file", plot_file);
	pp.query("plot_int", plot_int);
    }

    {
	ParmParse pp("adv");
	
	pp.query("cfl", cfl);
        pp.query("do_reflux", do_reflux);
    }
}

// set covered coarse cells to be the average of overlying fine cells
void
AmrAdv::AverageDown ()
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
	amrex::average_down(*phi_new[lev+1], *phi_new[lev],
			     geom[lev+1], geom[lev],
			     0, phi_new[lev]->nComp(), refRatio(lev));
    }
}

// more flexible version of AverageDown() that lets you average down across multiple levels
void
AmrAdv::AverageDownTo (int crse_lev)
{
    amrex::average_down(*phi_new[crse_lev+1], *phi_new[crse_lev],
			 geom[crse_lev+1], geom[crse_lev],
			 0, phi_new[crse_lev]->nComp(), refRatio(crse_lev));
}

// compute the number of cells at a level
long
AmrAdv::CountCells (int lev)
{
    const int N = grids[lev].size();

    long cnt = 0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:cnt)
#endif
    for (int i = 0; i < N; ++i)
    {
        cnt += grids[lev][i].numPts();
    }

    return cnt;
}

// compute a new multifab by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
AmrAdv::FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    if (lev == 0)
    {
	Array<MultiFab*> smf;
	Array<Real> stime;
	GetData(0, time, smf, stime);

	AmrAdvPhysBC physbc;
	amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
				     geom[lev], physbc);
    }
    else
    {
	Array<MultiFab*> cmf, fmf;
	Array<Real> ctime, ftime;
	GetData(lev-1, time, cmf, ctime);
	GetData(lev  , time, fmf, ftime);

	AmrAdvPhysBC cphysbc, fphysbc;
	Interpolater* mapper = &cell_cons_interp;

	int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR}; // periodic boundaryies
	int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
	Array<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

	amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
				   0, icomp, ncomp, geom[lev-1], geom[lev],
				   cphysbc, fphysbc, refRatio(lev-1),
				   mapper, bcs);
    }
}

// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void
AmrAdv::FillCoarsePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    BL_ASSERT(lev > 0);

    Array<MultiFab*> cmf;
    Array<Real> ctime;
    GetData(lev-1, time, cmf, ctime);
    
    if (cmf.size() != 1) {
	amrex::Abort("FillCoarsePatch: how did this happen?");
    }

    AmrAdvPhysBC cphysbc, fphysbc;
    Interpolater* mapper = &cell_cons_interp;
    
    int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR}; // periodic boundaryies
    int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    Array<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

    amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
				 cphysbc, fphysbc, refRatio(lev-1),
				 mapper, bcs);
}

// utility to copy in data from phi_old and/or phi_new into another multifab
void
AmrAdv::GetData (int lev, Real time, Array<MultiFab*>& data, Array<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
	data.push_back(phi_new[lev].get());
	datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
	data.push_back(phi_old[lev].get());
	datatime.push_back(t_old[lev]);
    }
    else
    {
	data.push_back(phi_old[lev].get());
	data.push_back(phi_new[lev].get());
	datatime.push_back(t_old[lev]);
	datatime.push_back(t_new[lev]);
    }
}

// in AmrAdvEvolve.cpp
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
            // update lev based on coarse-fine flux mismatch
	    flux_reg[lev+1]->Reflux(*phi_new[lev], 1.0, 0, 0, phi_new[lev]->nComp(), geom[lev]);
	}

	AverageDownTo(lev); // average lev+1 down to lev
    }
    
}

// advance a single level for a single time step, updates flux registers
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

            // compute velocities on faces (prescribed function of space and time)
	    get_face_velocity(lev, ctr_time,
			      AMREX_D_DECL(BL_TO_FORTRAN(uface[0]),
				     BL_TO_FORTRAN(uface[1]),
				     BL_TO_FORTRAN(uface[2])),
			      dx, prob_lo);

            // compute new state (stateout) and fluxes.
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
	    for (int i = 0; i < BL_SPACEDIM; ++i) {
		flux_reg[lev+1]->CrseInit(fluxes[i],i,0,0,fluxes[i].nComp(), -1.0);
	    }	    
	}
	if (flux_reg[lev]) {
	    for (int i = 0; i < BL_SPACEDIM; ++i) {
		flux_reg[lev]->FineAdd(fluxes[i],i,0,0,fluxes[i].nComp(), 1.0);
	    }
	}
    }
}

// a wrapper for EstTimeStep(0
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

// compute dt from CFL considerations
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

// get plotfile name
std::string
AmrAdv::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}

// put together an array of multifabs for writing
Array<const MultiFab*>
AmrAdv::PlotFileMF () const
{
    Array<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {
	r.push_back(phi_new[i].get());
    }
    return r;
}

// set plotfile variable names
Array<std::string>
AmrAdv::PlotFileVarNames () const
{
    return {"phi"};
}

// write plotfile to disk
void
AmrAdv::WritePlotFile () const
{
    const std::string& plotfilename = PlotFileName(istep[0]);
    const auto& mf = PlotFileMF();
    const auto& varnames = PlotFileVarNames();
    
    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
				    Geom(), t_new[0], istep, refRatio());
}
