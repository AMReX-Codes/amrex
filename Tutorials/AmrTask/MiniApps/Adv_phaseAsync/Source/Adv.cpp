#include <Adv.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>

#include <Perilla.H>
#include <RegionGraph.H>
#include <WorkerThread.H>

#include <iostream>
#include <iomanip>

#include <string>
#include <fstream>


using namespace amrex;
using namespace perilla;

int      Adv::verbose         = 0;
Real     Adv::cfl             = 0.9;
int      Adv::do_reflux       = 0;

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
}

Adv::Adv ()
{
    flux_reg = 0;
    Sborder = 0;
    SborderFPI.clear();    
    S_fine = NULL;
    S_crse = NULL;
    RG_S_fine = NULL;
    RG_S_crse = NULL;
}

Adv::Adv (Amr&            papa,
	int             lev,
	const Geometry& level_geom,
	const BoxArray& bl,
	const DistributionMapping& dm,
	Real            time)
    :
	AmrLevel(papa,lev,level_geom,bl,dm,time) 
{
    flux_reg = 0;
    if (level > 0 && do_reflux)
	flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);

    Sborder = 0;
    SborderFPI.clear();
    S_fine = NULL;
    S_crse = NULL;
    RG_S_fine = NULL;
    RG_S_crse = NULL;    
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

	initdata(level, cur_time, ARLIM_3D(lo), ARLIM_3D(hi),
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

    if(perilla::isMasterThread())
	if (do_reflux && level < finest_level)
	    reflux();
    perilla::syncAllWorkerThreads();

    avgDown(iteration);
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
    for (int k = finest_level; k>= 0; k--)
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

    Adv& fine_lev = getLevel(level+1);
    MultiFab&  S_fine   = fine_lev.get_new_data(State_Type);
    MultiFab&  S_crse   = get_new_data(State_Type);

    amrex::average_down(S_fine,S_crse,
	    fine_lev.geom,geom,
	    0,S_fine.nComp(),parent->refRatio(level));

}


    void
Adv::avgDown (int iteration)
{
    avgDown(State_Type, iteration);
}

    void
Adv::avgDown (int state_indx, int iteration)
{
    int tid = perilla::tid();
    int tg = perilla::wid();
    int nt = perilla::wtid();
    int myProc = ParallelDescriptor::MyProc();
    int f;

    if(level < parent->finestLevel())
    {
	Adv& fine_lev = getLevel(level+1);

	MultiFab* tS_fine = &(fine_lev.get_new_data(state_indx));
	MultiFab* tS_crse = &(get_new_data(state_indx));

	MultiFab& S_new = get_new_data(state_indx);
	const Real time = state[state_indx].curTime();


	if(perilla::isMasterWorkerThread())
	    RG_S_crse->Reset();

	while(!RG_S_crse->isGraphEmpty())
	{
	    f = RG_S_crse->getAnyFireableRegion();	  

	    Perilla::multifabCopyPull(RG_S_crse, RG_S_fine, tS_crse, tS_fine, f, 0, 0, S_fine->nComp(), 0, 0, false);

	    perilla::syncWorkerThreads();

	    // Send data to advance for next subcycle iteration  
	    if(iteration < parent->nCycle(level))
		SborderFPI[iteration]->FillPatchPush(NUM_GROW, time, state_indx, 0, NUM_STATE, f, 0x02, false);
	    RG_S_crse->finalizeRegion(f);
	}

	if(perilla::isMasterWorkerThread())
	    RG_S_crse->finalizeGraph();
    }

    perilla::syncWorkers();

    if (level > 0 && iteration == parent->nCycle(level))
    {
	Adv& crse_lev = getLevel(level-1);            
	MultiFab* tS_fine = &(get_new_data(state_indx));
	MultiFab* tS_crse = &(crse_lev.get_new_data(state_indx));


	if(perilla::isMasterWorkerThread())
	{
	    crse_lev.RG_S_fine->Reset();	  
	}


	while(!crse_lev.RG_S_fine->isGraphEmpty())
	{	
	    f = crse_lev.RG_S_fine->getAnyFireableRegion();
	    int lfi = tS_fine->IndexArray()[f];
#if 0
            const FArrayBox& fin = (*S_fine)[lfi]; 
            const FArrayBox& cin = (*crse_S_fine)[lfi]; 

	    int tg = perilla::wid();
	    int nt = perilla::wtid();

	    for(int t=0; t<crse_lev.RG_S_fine->fabTiles[f]->numTiles; t++)
		if(t % (perilla::NUM_THREADS_PER_TEAM-1) == nt)
		{
		    const Box& tbx = *(crse_lev.RG_S_fine->fabTiles[f]->tileBx[t]);
                    amrex_avgdown(tbx,cin.array(),fin.array(),0,0,ncomp,ratio);
		}
#endif
	    perilla::syncWorkerThreads();
	    Perilla::multifabCopyPushAsync(crse_lev.RG_S_crse, crse_lev.RG_S_fine, tS_crse, tS_fine, f, 0, 0, crse_lev.S_fine->nComp(), 0, 0, false);

	    crse_lev.RG_S_fine->finalizeRegion(f);	      
	}

	if(perilla::isMasterWorkerThread())
	    crse_lev.RG_S_fine->finalizeGraph();
    }
}

    void
Adv::initPerilla(Real time)
{
    int state_indx = State_Type;
    Sborder = new MultiFab(grids, dmap, NUM_STATE, NUM_GROW);

    SborderFPI.resize(parent->nCycle(level));        
    for(int i=0; i<parent->nCycle(level); i++)
    {
	SborderFPI[i] = new FillPatchIterator(*this, *(Sborder),NUM_GROW, time+(i*parent->dtLevel(level)), State_Type, 0, NUM_STATE,0);
	SborderFPI[i]->initFillPatch(NUM_GROW, time+(i*parent->dtLevel(level)), State_Type, 0, NUM_STATE, i+1);	
    }

    if(level < parent->finestLevel())
    {
	Adv& fine_lev = getLevel(level+1);
	S_fine = &(fine_lev.get_new_data(state_indx));
	S_crse = &(get_new_data(state_indx));

	const BoxArray& fine_BA = S_fine->boxArray();
	BoxArray crse_S_fine_BA = fine_BA;
	crse_S_fine_BA.coarsen(parent->refRatio(level));
	crse_S_fine = new MultiFab(crse_S_fine_BA, S_fine->DistributionMap(), S_fine->nComp(),0);

	RG_S_crse = new RegionGraph(S_crse->IndexArray().size());
	RG_S_fine = new RegionGraph(crse_S_fine->IndexArray().size());
	RG_S_fine->buildTileArray(*crse_S_fine);

	Perilla::multifabExtractCopyAssoc( RG_S_crse, RG_S_fine, *S_crse, *crse_S_fine, S_fine->nComp(), 0, 0, Periodicity::NonPeriodic());
	parent->graphArray[level].push_back(RG_S_crse);
	parent->graphArray[level].push_back(RG_S_fine);
    }
}

    void
Adv::finalizePerilla (Real time)
{
    if(ParallelDescriptor::MyProc()==0)
	std::cout<< "Finalizing Perilla Level " << level <<std::endl;

    for(int i=0; i< parent->nCycle(level); i++)
    {
	 delete SborderFPI[i];
    }
    SborderFPI.clear();
    if(Sborder) delete Sborder;    

    if(level < parent->finestLevel())
    {
	S_fine = 0;
	S_crse = 0;
	if(crse_S_fine) delete crse_S_fine;
	if(RG_S_fine) delete RG_S_fine;
	if(RG_S_crse) delete RG_S_crse;	
    }
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
	Vector<int>  itags;

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

	    state_error(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
		    BL_TO_FORTRAN_3D(S_new[mfi]),
		    &tagval, &clearval,
		    ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()),
		    ZFILL(dx), ZFILL(prob_lo), &time, &level);
	    // Now update the tags in the TagBox.
	    //
	    tagfab.tags_and_untags(itags, tilebx);
	}
    }
}
