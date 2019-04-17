#ifdef USE_PERILLA

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFabUtil_C.H>
#include <AMReX_MultiFabUtil_Perilla.H>

#include <Perilla.H>
#include <RegionGraph.H>
#include <WorkerThread.H>
#include <AMReX_Amr.H>
#include <AMReX_AmrLevel.H>

namespace amrex
{
  void average_down_push (RGIter& rgi, Amr& amr, MultiFab& S_fine, MultiFab& S_crse, MultiFab& crse_S_fine, RegionGraph* RG_fine, RegionGraph* RG_crse,
                       int scomp, int ncomp, const IntVect& ratio, int f)
    {
        if(rgi.currentItr != rgi.totalItr)
	  return;

	f = rgi.currentRegion;

        BL_ASSERT(S_crse.nComp() == S_fine.nComp());

	//if(RG_fine->graphID == 23)
	//std::cout<<"In avg down push f " << f << std::endl;

	//  NOTE: The tilebox is defined at the coarse level.
	int lfi = crse_S_fine.IndexArray()[f];
	//const Box& tbx = crse_S_fine[lfi].box();

	//  NOTE: We copy from component scomp of the fine fab into component 0 of the crse fab
	//        because the crse fab is a temporary which was made starting at comp 0, it is
	//        not part of the actual crse multifab which came in.

	int tg = perilla::wid();
	int nt = perilla::wtid();
	
	for(int t=0; t<RG_fine->fabTiles[f]->numTiles; t++)
	  if(t % (perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS) == nt-perilla::NUM_COMM_THREADS)
	    {
	      const Box& tbx = *(RG_fine->fabTiles[f]->tileBx[t]);
	      
              amrex_avgdown(tbx,crse_S_fine[lfi],S_fine[lfi],0,scomp,ncomp,ratio);
	    }
	RG_fine->worker[tg]->l_barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS); // Barrier to synchronize team threads
	
	//if(RG_fine->graphID == 23)
	//std::cout<<"In avg down pushAsych f " << f << std::endl;

	Perilla::multifabCopyPush( RG_crse, RG_fine, &S_crse, &crse_S_fine, f, scomp, 0, ncomp, 0, 0, false);
	
        //S_crse.copy(crse_S_fine,0,scomp,ncomp);
   }

  void average_down_pull (RGIter& rgi, MultiFab& S_fine, MultiFab& S_crse, RegionGraph* RG_fine, RegionGraph* RG_crse, 
                       int scomp, int ncomp, const IntVect& ratio, int f)
    {
        if(rgi.currentItr != 1)
	  return;

	f = rgi.currentRegion;

        BL_ASSERT(S_crse.nComp() == S_fine.nComp());

	Perilla::multifabCopyPull( RG_crse, RG_fine, &S_crse, &S_fine, f,scomp, 0, ncomp, 0, 0, false);
	
        //S_crse.copy(crse_S_fine,0,scomp,ncomp);
   }

  void average_down_push (RGIter& rgi, Amr& amr, MultiFab& S_fine, MultiFab& S_crse, MultiFab& crse_S_fine, RegionGraph* RG_fine, RegionGraph* RG_crse, 
			  const Geometry& fgeom, const Geometry& cgeom, 
			  int scomp, int ncomp, int rr, int f)
     {
       average_down_push(rgi,amr,S_fine,S_crse,crse_S_fine,RG_fine,RG_crse,fgeom,cgeom,scomp,ncomp,rr*IntVect::TheUnitVector(),f);
     }
  
  void average_down_pull (RGIter& rgi, MultiFab& S_fine, MultiFab& S_crse, RegionGraph* RG_fine, RegionGraph* RG_crse, const Geometry& fgeom, const Geometry& cgeom, 
		     int scomp, int ncomp, int rr, int f)
     {
       average_down_pull(rgi, S_fine,S_crse,RG_fine,RG_crse,fgeom,cgeom,scomp,ncomp,rr*IntVect::TheUnitVector(),f);
     }

  
  void average_down_push (RGIter& rgi, Amr& amr, MultiFab& S_fine, MultiFab& S_crse, MultiFab& crse_S_fine, RegionGraph* RG_fine, RegionGraph* RG_crse, 
			  const Geometry& fgeom, const Geometry& cgeom, 
			  int scomp, int ncomp, const IntVect& ratio, int f)
    {
        if(rgi.currentItr != rgi.totalItr)
	  return;

	f = rgi.currentRegion;

        if (S_fine.is_nodal() || S_crse.is_nodal())
        {
            amrex::Error("Can't use amrex::average_down for nodal MultiFab!");
        }

#if (BL_SPACEDIM == 3)
	amrex::average_down_push(rgi, amr, S_fine, S_crse, crse_S_fine, RG_fine, RG_crse, scomp, ncomp, ratio, f);
	return;
#else

        BL_ASSERT(S_crse.nComp() == S_fine.nComp());


	MultiFab fvolume;
	fgeom.GetVolume(fvolume, fine_BA, 0);

	//#ifdef _OPENMP
	//#pragma omp parallel
	//#endif
        //for (MFIter mfi(crse_S_fine,true); mfi.isValid(); ++mfi)
        //{
	//  NOTE: The tilebox is defined at the coarse level.
	int lfi = crse_S_fine.IndexArray()[f];
 	const Box& tbx = crse_S_fine[ lfi ].box();
	
        amrex_avgdown_with_vol(tbx,crse_S_fine[lfi],S_fine[lfi],fvolume[lfi],
                               0,scomp,ncomp,ratio);
	//}

	Perilla::multifabCopyPush( RG_crse, RG_fine, &S_crse, &crse_S_fine, f, scomp, 0, ncomp, 0, 0, false);
	
        //S_crse.copy(crse_S_fine,0,scomp,ncomp);
#endif
   }

  void average_down_pull (RGIter& rgi, MultiFab& S_fine, MultiFab& S_crse, RegionGraph* RG_fine, RegionGraph* RG_crse, const Geometry& fgeom, const Geometry& cgeom, 
			      int scomp, int ncomp, const IntVect& ratio, int f)
    {
        if(rgi.currentItr != 1)
	  return;

	f = rgi.currentRegion;

        if (S_fine.is_nodal() || S_crse.is_nodal())
        {
            amrex::Error("Can't use amrex::average_down for nodal MultiFab!");
        }

#if (BL_SPACEDIM == 3)
	amrex::average_down_pull(rgi, S_fine, S_crse, RG_fine, RG_crse, scomp, ncomp, ratio, f);
	return;
#else

        BL_ASSERT(S_crse.nComp() == S_fine.nComp());


	Perilla::multifabCopyPull(rgi, RG_crse, RG_fine, &S_crse, &S_fine, f, scomp, 0, ncomp, 0, 0, false);
	
        //S_crse.copy(crse_S_fine,0,scomp,ncomp);
#endif
   }
  
  void average_down_push (RGIter& rgi, Amr& amr, MultiFab& S_fine, MultiFab& S_crse, MultiFab& crse_S_fine, RegionGraph* RG_fine, RegionGraph* RG_crse,
			  int scomp, int ncomp, int rr, int f)
    {
      average_down_push(rgi,amr,S_fine,S_crse,crse_S_fine,RG_fine,RG_crse,scomp,ncomp,rr*IntVect::TheUnitVector(),f);
    }
  void average_down_pull (RGIter& rgi, MultiFab& S_fine, MultiFab& S_crse, RegionGraph* RG_fine, RegionGraph* RG_crse, int scomp, int ncomp, int rr, int f)
    {
      average_down_pull(rgi, S_fine,S_crse,RG_fine,RG_crse,scomp,ncomp,rr*IntVect::TheUnitVector(),f);
    }
}

#endif
