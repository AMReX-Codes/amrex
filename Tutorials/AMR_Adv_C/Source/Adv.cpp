
#include <Adv.H>
#include <Adv_F.H>
#include <VisMF.H>
#include <TagBox.H>
#include <ParmParse.H>

int      Adv::verbose         = 0;
Real     Adv::cfl             = 0.9;
int      Adv::do_reflux       = 1;

int      Adv::NUM_STATE       = 1;  // One variable in the state
int      Adv::NUM_GROW        = 3;  // number of ghost cells

// Note: Adv::variableSetUp is in Adv_setup.cpp

void
Adv::variableCleanUp () 
{
    desc_lst.clear();
}

void
Adv::read_params ()
{
    static bool done = false;

    if (done) return;

    done = true;

    ParmParse pp("adv");   

    pp.query("v",verbose);
    pp.query("cfl",cfl);
    pp.query("do_reflux",do_reflux);

    // This tutorial code only supports Cartesian coordinates.
    if (! Geometry::IsCartesian()) {
	BoxLib::Abort("Please set geom.coord_sys = 0");
    }

    // This tutorial code only supports periodic boundaries.
    if (! Geometry::isAllPeriodic()) {
	BoxLib::Abort("Please set geom.is_periodic = 1 1 1");
    }


}

Adv::Adv ()
{
    flux_reg = 0;
}

Adv::Adv (Amr&            papa,
	  int             lev,
	  const Geometry& level_geom,
	  const BoxArray& bl,
	  Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,time) 
{
    flux_reg = 0;
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,crse_ratio,level,NUM_STATE);
}

Adv::~Adv () 
{
    delete flux_reg;
}

void
Adv::initData ()
{
    //
    // Loop over grids, call FORTRAN function to init with data.
    //
    const Real* dx  = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    const Real* prob_hi = geom.ProbHi();
    MultiFab& S_new = get_new_data(State_Type);
    Real cur_time   = state[State_Type].curTime();

    if (verbose && ParallelDescriptor::IOProcessor())
	std::cout << "Initializing the data at level " << level << std::endl;

    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& box     = mfi.validbox();
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();

          BL_FORT_PROC_CALL(INITDATA,initdata)
	      (level, cur_time, ARLIM_3D(lo), ARLIM_3D(hi),
	       BL_TO_FORTRAN_3D(S_new[mfi]), ZFILL(dx),
	       ZFILL(prob_lo), ZFILL(prob_hi));
    }

    if (verbose && ParallelDescriptor::IOProcessor())
	std::cout << "Done initializing the level " << level << " data " << std::endl;
}

void
Adv::init (AmrLevel &old)
{
    Adv* oldlev = (Adv*) &old;
    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev->state[State_Type].curTime();
    Real prev_time = oldlev->state[State_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(State_Type);
    
    for (FillPatchIterator fpi(old,S_new,0,cur_time,State_Type,0,NUM_STATE);
          fpi.isValid();
          ++fpi)
    {
        S_new[fpi].copy(fpi());
    }
}

//
// This version inits the data on a new level that did not
// exist before regridding.
//
void
Adv::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[State_Type].curTime();
    Real prev_time = getLevel(level-1).state[State_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);
    MultiFab& S_new = get_new_data(State_Type);
    FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, NUM_STATE);
}

Real
Adv::initialTimeStep ()
{
    return estTimeStep(0.0);
}

Real
Adv::estTimeStep (Real)
{
    // This is just a dummy value to start with 
    Real dt_est  = 1.0e+20;

    const Real* dx = geom.CellSize();
    const MultiFab& stateMF = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_est)
#endif
    {
	for (MFIter mfi(stateMF, true); mfi.isValid(); ++mfi)
	{
	    const Box& box = mfi.tilebox();

	    BL_FORT_PROC_CALL(ESTDT,estdt)
		(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
		 BL_TO_FORTRAN_3D(stateMF[mfi]),
		 ZFILL(dx), &dt_est);
	}
    }

    ParallelDescriptor::ReduceRealMin(dt_est);
    dt_est *= cfl;

    if (verbose && ParallelDescriptor::IOProcessor())
	std::cout << "Adv::estTimeStep at level " << level << ":  dt_est = " << dt_est << std::endl;
    
    return dt_est;
}

void
Adv::computeNewDt (int                   finest_level,
		   int                   sub_cycle,
		   Array<int>&           n_cycle,
		   const Array<IntVect>& ref_ratio,
		   Array<Real>&          dt_min,
		   Array<Real>&          dt_level,
		   Real                  stop_time,
		   int                   post_regrid_flag)
{
    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;

    for (int i = 0; i <= finest_level; i++)
    {
        Adv& adv_level = getLevel(i);
        dt_min[i] = adv_level.estTimeStep(dt_level[i]);
    }

    if (post_regrid_flag == 1) 
    {
	//
	// Limit dt's by pre-regrid dt
	//
	for (int i = 0; i <= finest_level; i++)
	{
	    dt_min[i] = std::min(dt_min[i],dt_level[i]);
	}
    } 
    
    //
    // Find the minimum over all levels
    //
    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

void
Adv::computeInitialDt (int                   finest_level,
		       int                   sub_cycle,
		       Array<int>&           n_cycle,
		       const Array<IntVect>& ref_ratio,
		       Array<Real>&          dt_level,
		       Real                  stop_time)
{
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor   *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

void
Adv::post_timestep (int iteration)
{
    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level)
        reflux();

    if (level < finest_level)
        avgDown();
}

void
Adv::post_init (Real stop_time)
{
    if (level > 0)
        return;
    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level-1; k>= 0; k--)
        getLevel(k).avgDown();
}

void
Adv::reflux ()
{
    BL_ASSERT(level<parent->finestLevel());

    const Real strt = ParallelDescriptor::second();

    getFluxReg(level+1).Reflux(get_new_data(State_Type),1.0,0,0,NUM_STATE,geom);
    
    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;
	
        ParallelDescriptor::ReduceRealMax(end,IOProc);
	
        if (ParallelDescriptor::IOProcessor())
            std::cout << "Adv::reflux() at level " << level << " : time = " << end << std::endl;
    }
}

void
Adv::avgDown ()
{
    if (level == parent->finestLevel()) return;
    avgDown(State_Type);
}

void
Adv::avgDown (int state_indx)
{
    if (level == parent->finestLevel()) return;

    Adv& fine_lev = getLevel(level+1);
    MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);
    MultiFab&  S_crse   = get_new_data(state_indx);
    
    BoxLib::average_down(S_fine,S_crse,
                         fine_lev.geom,geom,
                         0,S_fine.nComp(),parent->refRatio(level));
}

void
Adv::errorEst (TagBoxArray& tags,
	       int          clearval,
	       int          tagval,
	       Real         time,
	       int          n_error_buf,
	       int          ngrow)
{
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<int>  itags;
	
	for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
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

	    BL_FORT_PROC_CALL(STATE_ERROR, state_error)
		(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
		 BL_TO_FORTRAN_3D(S_new[mfi]),
		 &tagval, &clearval, 
		 ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()), 
		 ZFILL(dx), ZFILL(prob_lo), &time, &level);
	}
    }
}
