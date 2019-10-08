#include <Adv.H>
#include <Adv_F.H>
#include <WorkerThread.H>
#include <RegionGraph.H>
#include <Perilla.H>
#include <string>
#include <fstream>

using namespace perilla;

Real
Adv::advance (Real time,
	Real dt,
	int  iteration,
	int  ncycle)
{
    int tg = perilla::wid();
    int nt = perilla::wtid();
    int myProc = ParallelDescriptor::MyProc();

    perilla::syncAllWorkerThreads(); //tasks handled by workers in a process share the state data
    if(perilla::isMasterThread())
    {
	for (int k = 0; k < NUM_STATE_TYPE; k++) {
	    state[k].allocOldData();
	    state[k].swapTimeLevels(dt);
	}
    }
    perilla::syncAllWorkerThreads();

    MultiFab& S_new = get_new_data(State_Type);

    const Real prev_time = state[State_Type].prevTime();
    const Real cur_time = state[State_Type].curTime();
    const Real ctr_time = 0.5*(prev_time + cur_time);

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    // Get pointers to Flux registers, or set pointer to zero if not there.
    FluxRegister *fine    = 0;
    FluxRegister *current = 0;

    int finest_level = parent->finestLevel();

    if(perilla::isMasterThread())
    {
	if (do_reflux && level < finest_level) {
	    fine = &getFluxReg(level+1);
	    fine->setVal(0.0);
	}
    }
    perilla::syncAllWorkerThreads();

    if (do_reflux && level > 0) {
	current = &getFluxReg(level);
    }

    MultiFab fluxes[BL_SPACEDIM];

    if (do_reflux)
    {
	for (int j = 0; j < BL_SPACEDIM; j++)
	{
	    BoxArray ba = S_new.boxArray();
	    ba.surroundingNodes(j);
	    fluxes[j].define(ba, dmap, NUM_STATE, 0);
	}
    }

    //This is the FillPatch on the coarsest level, starting the async execution on multiple AMR levels
    if(level == 0 && iteration == 1)
    {
	Adv& thisLevel = getLevel(level);	
	bool isSingleThread=false;
	for(int f=0; f<Sborder->IndexArray().size(); f++)
	{
	    if(WorkerThread::isMyRegion(tg, f))
	    {		
		for(int i=0; i<parent->nCycle(level); i++)
		{
		    //cout<<"Filling level 0"<<endl;
		    thisLevel.SborderFPI[i]->FillPatchPush(NUM_GROW, time+(i*parent->dtLevel(level)), State_Type /*index to state data*/, 0 /*base component*/, NUM_STATE/* num components*/, f /*fill FAB f only*/, 0xFF, isSingleThread);				
		}
	    }
	}
    }    
    perilla::syncWorkerThreads();

    FArrayBox flux[BL_SPACEDIM], uface[BL_SPACEDIM];

    if(level < parent->finestLevel())
    {
	Adv& upperLevel = getLevel(level+1);
	bool isSingleThread=false;
	unsigned char tuc = 0x04;
	upperLevel.SborderFPI[0]->FillPatchPush(NUM_GROW, time, State_Type /*index to state data*/, 0 /*base component*/, NUM_STATE/* num components*/, -1/*fill whole AMR level*/, tuc, isSingleThread);
    }

    //double ulp_etime = omp_get_wtime();

    if(perilla::isMasterWorkerThread())
	SborderFPI[iteration-1]->Reset();
    perilla::syncWorkerThreads();

    double stime, etime;

    if(nt == perilla::NUM_THREADS_PER_TEAM-2)
    {
	bool singleThread= true;
	while(SborderFPI[iteration-1]->destGraph->worker[tg]->completedRegionQueue->queueSize(true)!=SborderFPI[iteration-1]->destGraph->worker[tg]->totalTasks ||
		SborderFPI[iteration-1]->destGraph->worker[tg]->computedTasks!=SborderFPI[iteration-1]->destGraph->worker[tg]->totalTasks)
	{
	    int f = SborderFPI[iteration-1]->destGraph->getFireableRegion(singleThread);
	    if(f != -1)
	    {
		SborderFPI[iteration-1]->FillPatchPull(NUM_GROW, time, State_Type, 0, NUM_STATE, f, true);
		SborderFPI[iteration-1]->destGraph->setFireableRegion(f);
	    }

	    if(SborderFPI[iteration-1]->destGraph->worker[tg]->computedRegionQueue->queueSize()!=0)
	    {
		f = SborderFPI[iteration-1]->destGraph->worker[tg]->computedRegionQueue->removeRegion();

		if(level == parent->finestLevel() && iteration < ncycle)
		    SborderFPI[iteration]->FillPatchPush(NUM_GROW, time+dt, State_Type, 0, NUM_STATE, f /*fill FAB f only*/, 0x02, true);

		if(level < parent->finestLevel())
		{
		    Adv& upperLevel = getLevel(level+1);
		    for(int i=0; i<parent->nCycle(level+1); i++)
		    {
			unsigned char tuc = 0x01;
			upperLevel.SborderFPI[i]->FillPatchPush(NUM_GROW, time+(i*parent->dtLevel(level+1)), State_Type, 0, NUM_STATE, f /*fill FAB f only*/, tuc, true);
		    }
		}	
		SborderFPI[iteration-1]->destGraph->worker[tg]->completedRegionQueue->addRegion(f,true);
	    }
	}
    }
    else
    {
	while(!SborderFPI[iteration-1]->destGraph->isGraphEmptyV2())
	{
	    int f = SborderFPI[iteration-1]->destGraph->getPulledFireableRegion();	
	    int fid = S_new.IndexArray()[f];
	    int fis = Sborder->IndexArray()[f];
	    const FArrayBox& statein = (*(Sborder))[fis];
	    FArrayBox& stateout      =   S_new[fid];
	    MFIter mfi(S_new, true);
            SborderFPI[iteration-1]->destGraph->worker[tg]->l_barr->sync(perilla::NUM_THREADS_PER_TEAM-2);
	    for(int t=0; t<SborderFPI[iteration-1]->destGraph->fabTiles[f]->numTiles; t++)	  
		if(t % (perilla::NUM_THREADS_PER_TEAM-2) == nt)
		{	
		    const Box& bx = *(SborderFPI[iteration-1]->destGraph->fabTiles[f]->tileBx[t]);
		    for (int i = 0; i < BL_SPACEDIM ; i++) {
			const Box& bxtmp = amrex::surroundingNodes(bx,i);
			{
			    flux[i].resize(bxtmp,NUM_STATE);
			    uface[i].resize(amrex::grow(bxtmp,1),1);
			}
		    }

		    get_face_velocity(level, ctr_time,
			    D_DECL(BL_TO_FORTRAN(uface[0]),
				BL_TO_FORTRAN(uface[1]),
				BL_TO_FORTRAN(uface[2])),
			    dx, prob_lo);
		    advect(time, bx.loVect(), bx.hiVect(),
			    BL_TO_FORTRAN_3D(statein),
			    BL_TO_FORTRAN_3D(stateout),
			    D_DECL(BL_TO_FORTRAN_3D(uface[0]),
				BL_TO_FORTRAN_3D(uface[1]),
				BL_TO_FORTRAN_3D(uface[2])),
			    D_DECL(BL_TO_FORTRAN_3D(flux[0]),
				BL_TO_FORTRAN_3D(flux[1]),
				BL_TO_FORTRAN_3D(flux[2])),
			    dx, dt);

		    if (do_reflux) {
			for (int i = 0; i < BL_SPACEDIM ; i++)
			    fluxes[i][f].copy(flux[i],mfi.nodaltilebox(i));	  
		    }
		}   // End of single Fab computation using tiling
            SborderFPI[iteration-1]->destGraph->worker[tg]->l_barr->sync(perilla::NUM_THREADS_PER_TEAM-2);
	    if (do_reflux) {
		if (current) {
		    for (int i = 0; i < BL_SPACEDIM ; i++)
			current->FineAdd(fluxes[i][f], i, fluxes[i].IndexArray()[f], 0, 0, NUM_STATE, 1.0, RunOn::Cpu);
		}
	    }
	    SborderFPI[iteration-1]->destGraph->regionComputed(f);
	}// while(!SborderFPI[i]->destGraph->isGraphEmpty())
    }
    if (do_reflux) {
	if (fine) {
	    if(perilla::isMasterThread())
		for (int i = 0; i < BL_SPACEDIM ; i++)
		    fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.);
	}
	perilla::syncAllWorkerThreads();
    }

    if(perilla::isMasterWorkerThread())
	SborderFPI[iteration-1]->finalizeGraphs();

    perilla::syncAllWorkerThreads();
    return dt;
}
