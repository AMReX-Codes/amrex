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

    amrex_avgdown_with_vol(tbx,crse_S_fine[lfi],S_fine[lfi],fvolume[mfi],
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
    RG_fine->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads
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

    f = rgi.currentRegion;
    //  NOTE: The tilebox is defined at the coarse level.
    int lfi = crse_S_fine->IndexArray()[f];

    //  NOTE: We copy from component scomp of the fine fab into component 0 of the crse fab
    //        because the crse fab is a temporary which was made starting at comp 0, it is
    //        not part of the actual crse multifab which came in.

    perilla::syncWorkerThreads();
    int nThreads= perilla::nWorkerThreads();
    for(int t=0; t<RG_fine->fabTiles[f]->numTiles; t+= nThreads)
    {
        const Box& tbx = *(RG_fine->fabTiles[f]->tileBx[t]);
        amrex_avgdown(tbx,(*crse_S_fine)[lfi],(*S_fine)[lfi],0,scomp,ncomp,ratio);
    }
    perilla::syncWorkerThreads();
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


// *************************************************************************************************************

#if 0
// Average fine face-based MultiFab onto crse fine-centered MultiFab.
// This routine assumes that the crse BoxArray is a coarsened version of the fine BoxArray.
void average_down_faces (PArray<MultiFab>& fine, PArray<MultiFab>& crse, IntVect& ratio)
{
    BL_ASSERT(crse.size()  == BL_SPACEDIM);
    BL_ASSERT(fine.size()  == BL_SPACEDIM);
    BL_ASSERT(crse[0].nComp() == fine[0].nComp());

    int ncomp = crse[0].nComp();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (int n=0; n<BL_SPACEDIM; ++n) {
	for (MFIter mfi(crse[n],true); mfi.isValid(); ++mfi)
	{
	    const Box& tbx = mfi.tilebox();

	    BL_FORT_PROC_CALL(BL_AVGDOWN_FACES,bl_avgdown_faces)
		(tbx.loVect(),tbx.hiVect(),
		 BL_TO_FORTRAN(fine[n][mfi]),
		 BL_TO_FORTRAN(crse[n][mfi]),
		 ratio.getVect(),n,ncomp);
	}
    }
}
#endif
