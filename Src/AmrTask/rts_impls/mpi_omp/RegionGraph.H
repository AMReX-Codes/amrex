#ifndef P_REGIONGRAPH_H
#define P_REGIONGRAPH_H

#include <RegionQueue.H>
#include <LocalConnection.H>
#include <RemoteConnection.H>
#include <WorkerThread.H>
#include <AMReX_MultiFab.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>
#include <omp.h>

using namespace perilla;
#ifdef PERILLA_DEBUG
#include "PerillaMemCheck.H"
extern PerillaMemCheck memcheck;
#endif

namespace amrex{

    class RegionGraph;

    class FabCopyAssoc
    {
	public:
	    LocalConnection l_con;
	    RemoteConnection r_con;
	    FabCopyAssoc *next;
	    FabCopyAssoc *prev;
	    RegionGraph *graphPartner;
	    //  CopyAssoc *orig_copyAssoc;
	    FabCopyAssoc()
	    {
		next=0;
		prev=0;
		graphPartner=0;
	    }
            ~FabCopyAssoc()
            {
                if(next) delete next;
            }
    };

    class pTileArray
    {
	public:
	    int numTiles;
	    std::vector<Box*> tileBx;
	    Box* validBx;
	    pTileArray(): numTiles(0), tileBx(0){}
	    ~pTileArray()
	    {
		tileBx.clear();
	    }
    };

    class Task
    {
	public:
            std::vector<int> state;
	    bool init;
	    FabCopyAssoc *cpAsc_srcHead;
	    FabCopyAssoc *cpAsc_dstHead;

	    std::vector<int> depTaskIDs;
	    int numDepTasksCompleted;
	    bool depTasksCompleted;

	    Task()
	    {
		state.reserve(16);
		depTaskIDs.reserve(1024);
		depTasksCompleted = true;
		numDepTasksCompleted = 0;
		cpAsc_srcHead=0;
		cpAsc_dstHead=0;
	    }

	    ~Task()
	    {
		state.clear();
		depTaskIDs.clear();
		if(cpAsc_srcHead != 0)
		    delete[] cpAsc_srcHead;
		if(cpAsc_dstHead != 0)
		    delete[] cpAsc_dstHead;
	    }
    };

    class Worker
    {
	public:
	    int totalTasks;
	    int computedTasks;
	    bool init;
	    Barrier* barr;
	    Barrier* l_barr;
	    RegionQueue* fireableRegionQueue;
	    RegionQueue* unfireableRegionQueue;
	    RegionQueue* computedRegionQueue;
	    RegionQueue* completedRegionQueue;
	    Worker():init(false), l_barr(0), barr(0), totalTasks(0){}

            ~Worker(){
                delete barr;
                delete l_barr;
                delete fireableRegionQueue;
                delete unfireableRegionQueue;
                delete computedRegionQueue;
                delete completedRegionQueue;
            }
    };

    //template <class FAB>
    class CopyMap
    {
	public:
	    std::vector<FArrayBox*> map;
	    CopyMap *next;
	    CopyMap()
		:
		    next(0)
	{
#ifdef PERILLA_DEBUG
                memcheck.add(memcheck.genKey(this), (void*)this, "CopyMap");
#endif

        }

	    void alloc_CopyMap(const MultiFab& mf)
	    {

		Vector<int> IndArr = mf.IndexArray();
		const int n = IndArr.size();
		//const int n = mf.size();
		map.reserve(n);
		//sMap.resize(n);
		for (int i = 0; i < n; ++i)
		{
		    int K = IndArr[i];
		    const Box& tmp = mf.fabbox(K);
		    map.push_back(new FArrayBox(tmp, mf.nComp(), false, true));
		}
	    }
	    ~CopyMap()
	    {
                for (int i = 0; i < map.size(); ++i)
                {
                    delete map[i];
                }
		map.clear();
		if(next !=0 )
		    delete next;
#ifdef PERILLA_DEBUG
                memcheck.remove(memcheck.genKey(this));
#endif
	    }
    };

    class RegionGraph
    {
	public:
	    static int graphCnt;
	    int graphID;
	    int numTasks;
	    int numFabs;
	    int totalFinishes;
	    bool isDepGraph;
	    bool* okToReset;
	    omp_lock_t finishLock;

	    std::vector<pTileArray*> fabTiles;
	    std::vector<pTileArray*> fabTiles_gtbx;

	    std::vector<FArrayBox*> lMap;
	    std::vector<FArrayBox*> sMap;
	    std::vector<FArrayBox*> rMap;
	    CopyMap *sCopyMapHead;
	    CopyMap *rCopyMapHead;  

	    std::vector<Task*> task;
	    std::vector<Worker*> worker;

	    RegionGraph* srcLinkGraph;

	public:
	    RegionGraph(int numtasks);
	    void Initialize();
	    void Reset();
	    bool isGraphEmpty();
	    bool isGraphEmptyV2();
	    void finalizeGraph();
	    void regionGraphReset(int numfabs);
	    void regionGraphMinReset(void);
	    void enableAllRegions();
	    void disableRegion(int r, int tg);
	    void finalizeRegion(int r);
            void finalizeRegionGraph();
	    void regionComputed(int r);
	    bool isFireableRegion(int r);
	    int  getAnyFireableRegion();
            int  getAnyFireableRegion(RegionGraph& depGraph);
	    int  getPulledFireableRegion();
	    int  getFireableRegion(bool isSingleThread=false);
	    void setFireableRegion(int r);
	    void graphTeardown();
	    void workerTeardown();
            int size(){return task.size();}

	    int getRegion(){
	        return worker[perilla::wid()]->computedRegionQueue->getFrontRegion(true);
            }

            void syncComputeWorkerThreads(){
	        worker[perilla::wid()]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
            }
            void syncComputeWorkerThreads(int numthreads){
	        worker[perilla::wid()]->barr->sync(numthreads);
            }

            void syncWorkerThreads(){
	        worker[perilla::wid()]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
            }
            void syncWorkerThreads(int numthreads){
	        worker[perilla::wid()]->barr->sync(numthreads);
            }

	    void alloc_lMap(const MultiFab& mf)
	    {
		const int n = mf.IndexArray().size();
		lMap.reserve(n);
		for (int i = 0; i < n; ++i)
		{
		    int K = mf.IndexArray()[i];
		    const Box& tmp = mf.fabbox(K);
		    lMap.push_back(new FArrayBox(tmp, mf.nComp(), false, true));
		}
	    }  

	    void alloc_sMap(const MultiFab& mf)
	    {
		const int n = mf.IndexArray().size();
		sMap.reserve(n);
		for (int i = 0; i < n; ++i)
		{
		    int K = mf.IndexArray()[i];
		    const Box& tmp = mf.fabbox(K);
		    sMap.push_back(new FArrayBox(tmp, mf.nComp(), false, true));
		}
	    }

	    void alloc_rMap(const MultiFab& mf)
	    {
		const int n = mf.IndexArray().size();
		rMap.reserve(n);
		for (int i = 0; i < n; ++i)
		{
		    int K = mf.IndexArray()[i];
		    const Box& tmp = mf.fabbox(K);
		    rMap.push_back(new FArrayBox(tmp, mf.nComp(), false, true));
		}
	    }

	    void buildTileArray(const MultiFab& mf)
	    {
		const int n = mf.IndexArray().size();
		fabTiles.resize(n);

		for (int i = 0; i < n; ++i)
		{
		    fabTiles[i] = new pTileArray();
		}
		for (MFIter mfi(mf, true); mfi.isValid(); ++mfi)
		{
		    fabTiles[mfi.LocalIndex()]->numTiles++;
		    fabTiles[mfi.LocalIndex()]->tileBx.push_back(new Box(mfi.tilebox()));
		    fabTiles[mfi.LocalIndex()]->validBx = new Box(mfi.validbox());
		}
	    }

	    void buildTileArray_gtbx(const MultiFab& mf, int ng)
	    {
		const int n = mf.IndexArray().size();
		fabTiles_gtbx.resize(n);

		for (int i = 0; i < n; ++i)
		{
		    fabTiles_gtbx[i] = new pTileArray();
		}
		for (MFIter mfi(mf, true); mfi.isValid(); ++mfi)
		{
		    fabTiles_gtbx[mfi.LocalIndex()]->numTiles++;
		    fabTiles_gtbx[mfi.LocalIndex()]->tileBx.push_back(new Box(mfi.growntilebox(ng)));
		}
	    }

            void buildTileArray(const amrex::MultiFab& mf, const amrex::IntVect& tilesize)
            {
                int myProc = amrex::ParallelDescriptor::MyProc();
                const int n = mf.indexArray.size();
                fabTiles.resize(n);

                //typ = mf.boxArray().ixType();

                for (int i = 0; i < n; ++i)
                {
                    fabTiles[i] = new pTileArray();
                }
                for (amrex::MFIter mfi(mf, tilesize); mfi.isValid(); ++mfi)
                {
                    if( fabTiles[mfi.LocalIndex()]->numTiles == 0 )
                        fabTiles[mfi.LocalIndex()]->validBx = new amrex::Box(mfi.validbox());
                    fabTiles[mfi.LocalIndex()]->numTiles++;
                    fabTiles[mfi.LocalIndex()]->tileBx.push_back(new amrex::Box(mfi.tilebox()));
                }
            }


	    ~RegionGraph();
    };
}//end namespace


#endif
