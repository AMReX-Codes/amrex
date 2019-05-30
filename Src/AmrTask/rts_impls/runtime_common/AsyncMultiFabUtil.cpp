#include <AMReX_MultiFabUtil.H>
//#include <AMReX_MultiFabUtil_F.H>
#include <AsyncMultiFabUtil.H>
#include <AMReX_MultiFabUtil_C.H>
#include <Perilla.H>
#include <WorkerThread.H>

using namespace amrex;
using namespace perilla;

void average_down_push (Amr& amr, MultiFab& S_fine, MultiFab& S_crse, MultiFab& crse_S_fine, RegionGraph* RG_fine, RegionGraph* RG_crse, 
	const Geometry& fgeom, const Geometry& cgeom, int scomp, int ncomp, int rr, int f)
{
    average_down_push(amr,S_fine,S_crse,crse_S_fine,RG_fine,RG_crse,fgeom,cgeom,scomp,ncomp,rr*IntVect::TheUnitVector(),f);
}

void average_down_pull (MultiFab& S_fine, MultiFab& S_crse, RegionGraph* RG_fine, RegionGraph* RG_crse, 
	const Geometry& fgeom, const Geometry& cgeom, int scomp, int ncomp, int rr, int f)
{
    average_down_pull(S_fine,S_crse,RG_fine,RG_crse,fgeom,cgeom,scomp,ncomp,rr*IntVect::TheUnitVector(),f);
}

void average_down_push (Amr& amr, MultiFab& S_fine, MultiFab& S_crse, MultiFab& crse_S_fine, RegionGraph* RG_fine, RegionGraph* RG_crse, 
	const Geometry& fgeom, const Geometry& cgeom, int scomp, int ncomp, const IntVect& ratio, int f)
{
    if (S_fine.is_nodal() || S_crse.is_nodal())
    {
	amrex::Error("Can't use amrex::average_down for nodal MultiFab!");
    }

#if (BL_SPACEDIM == 3)
    average_down_push(amr, S_fine, S_crse, crse_S_fine, RG_fine, RG_crse, scomp, ncomp, ratio, f);
    return;
#else

    assert(S_crse.nComp() == S_fine.nComp());


    MultiFab fvolume;
    fgeom.GetVolume(fvolume, fine_BA, 0);

    int lfi = crse_S_fine.IndexArray()[f];
    const Box& tbx = crse_S_fine[ lfi ].box();

    amrex_avgdown_with_vol(tbx,crse_S_fine[lfi],S_fine[lfi],fvolume[lfi],
	    0,scomp,ncomp,ratio);

    Perilla::multifabCopyPushAsync(RG_crse, RG_fine, &S_crse, &crse_S_fine, f, scomp, 0, ncomp, 0, 0, false);
#endif
}

void average_down_pull (MultiFab& S_fine, MultiFab& S_crse, RegionGraph* RG_fine, RegionGraph* RG_crse, const Geometry& fgeom, const Geometry& cgeom, 
	int scomp, int ncomp, const IntVect& ratio, int f)
{

    if (S_fine.is_nodal() || S_crse.is_nodal())
    {
	amrex::Error("Can't use amrex::average_down for nodal MultiFab!");
    }

#if (BL_SPACEDIM == 3)
    average_down_pull(S_fine, S_crse, RG_fine, RG_crse, scomp, ncomp, ratio, f);
    return;
#else
    assert(S_crse.nComp() == S_fine.nComp());
    Perilla::multifabCopyPull(RG_crse, RG_fine, &S_crse, &S_fine, f, scomp, 0, ncomp, 0, 0, false);
#endif
}

// *************************************************************************************************************

void average_down_push (Amr& amr, MultiFab& S_fine, MultiFab& S_crse, MultiFab& crse_S_fine, RegionGraph* RG_fine, RegionGraph* RG_crse,
	int scomp, int ncomp, int rr, int f)
{
    average_down_push(amr,S_fine,S_crse,crse_S_fine,RG_fine,RG_crse,scomp,ncomp,rr*IntVect::TheUnitVector(),f);
}

void average_down_pull (MultiFab& S_fine, MultiFab& S_crse, RegionGraph* RG_fine, RegionGraph* RG_crse, int scomp, int ncomp, int rr, int f)
{
    average_down_pull(S_fine,S_crse,RG_fine,RG_crse,scomp,ncomp,rr*IntVect::TheUnitVector(),f);
}

void average_down_push (Amr& amr, MultiFab& S_fine, MultiFab& S_crse, MultiFab& crse_S_fine, RegionGraph* RG_fine, RegionGraph* RG_crse,
	int scomp, int ncomp, const IntVect& ratio, int f)
{
    assert(S_crse.nComp() == S_fine.nComp());

    //  NOTE: The tilebox is defined at the coarse level.
    int lfi = crse_S_fine.IndexArray()[f];
    int tg = WorkerThread::perilla_wid();
    int nt = WorkerThread::perilla_wtid();

    for(int t=0; t<RG_fine->fabTiles[f]->numTiles; t++)
	if(t % (perilla::NUM_THREADS_PER_TEAM-1) == nt)
	{
	    const Box& tbx = *(RG_fine->fabTiles[f]->tileBx[t]);
	    amrex_avgdown(tbx,crse_S_fine[lfi],S_fine[lfi],0,scomp,ncomp,ratio);
	}
    RG_fine->worker[tg]->l_barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads
    Perilla::multifabCopyPushAsync(RG_crse, RG_fine, &S_crse, &crse_S_fine, f, scomp, 0, ncomp, 0, 0, false);
}

void average_down_pull (MultiFab& S_fine, MultiFab& S_crse, RegionGraph* RG_fine, RegionGraph* RG_crse, 
	int scomp, int ncomp, const IntVect& ratio, int f)
{
    assert(S_crse.nComp() == S_fine.nComp());
    Perilla::multifabCopyPull(RG_crse, RG_fine, &S_crse, &S_fine, f, scomp, 0, ncomp, 0, 0, false);
}


void average_down_push (RGIter& rgi, MultiFab* S_fine, MultiFab* S_crse, MultiFab* crse_S_fine, RegionGraph* RG_fine, RegionGraph* RG_crse,amrex::Geometry& geom, amrex::Geometry& geom1,
	int scomp, int ncomp, const IntVect& ratio, int f)
{
    if(rgi.currentItr != rgi.totalItr)
	return;
    int tg = WorkerThread::perilla_wid();

    f = rgi.currentRegion;
    //  NOTE: The tilebox is defined at the coarse level.
    int lfi = crse_S_fine->IndexArray()[f];

    //  NOTE: We copy from component scomp of the fine fab into component 0 of the crse fab
    //        because the crse fab is a temporary which was made starting at comp 0, it is
    //        not part of the actual crse multifab which came in.

    //perilla::syncWorkerThreads();
    RG_fine->worker[tg]->l_barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads
    int nThreads= perilla::nWorkerThreads();
    for(int t=0; t<RG_fine->fabTiles[f]->numTiles; t+= nThreads)
    {
	const Box& tbx = *(RG_fine->fabTiles[f]->tileBx[t]);
	amrex_avgdown(tbx,(*crse_S_fine)[lfi],(*S_fine)[lfi],0,scomp,ncomp,ratio);
    }
    //perilla::syncWorkerThreads();
    RG_fine->worker[tg]->l_barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads
    Perilla::multifabCopyPush(RG_crse, RG_fine, S_crse, crse_S_fine, f, scomp, 0, ncomp, 0, 0, false);
}

void average_down_pull (RGIter& rgi, MultiFab* S_fine, MultiFab* S_crse, RegionGraph* RG_fine, RegionGraph* RG_crse, amrex::Geometry& geom, amrex::Geometry& geom1,
	int scomp, int ncomp, const IntVect& ratio, int f)
{
    if(rgi.currentItr != 1)
	return;
    f = rgi.currentRegion;

    Perilla::multifabCopyPull(RG_crse, RG_fine, S_crse, S_fine, f, scomp, 0, ncomp, 0, 0, false);
}


#if 0
#include "PerillaMemCheck.H"

void PerillaMemCheck::add(string key, void* obj, string classname)
{
    lock.lock();
    if(objMap.find(key) == objMap.end())
    {
        objMap[key]= obj;
	printf("Adding an object\n");
    }
    else{
        printf("Reinsert an object\n");
        exit(0);
    }
    lock.unlock();
}


void PerillaMemCheck::remove(string key){
    lock.lock();
    if(objMap.find(key) != objMap.end())
    {
        objMap.erase(key);
	printf("Removing an object\n");
    }
    else{
        printf("Object not found\n");
        exit(0);
    }

    lock.unlock();
}
void PerillaMemCheck::report(){
    if(objMap.size()) {
        printf("Memory leak found\n");
    }else printf("all packages deallocated\n");
}


#endif
