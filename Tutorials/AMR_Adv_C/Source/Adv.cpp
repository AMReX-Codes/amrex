
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
	       ZFILL(prob_lo));
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

    FillPatch(old, S_new, 0, cur_time, State_Type, 0, NUM_STATE);
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
	    //
	    // Now update the tags in the TagBox.
	    //
	    tagfab.tags_and_untags(itags, tilebx);
	}
    }
}
