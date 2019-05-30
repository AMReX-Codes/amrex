#ifndef _PERILLA_
#define _PERILLA_

#include <AMReX_MultiFab.H>
#include <vector>
#include <map>
#include <RegionGraph.H>
#include <pthread.h>

#include <RGIter.H>

//#define USE_PERILLA_PTHREADS

using namespace std;

namespace amrex{
    class Perilla
    {
	static int tagGen(int src, int dest, int channelID, int nFabs, int nChannels);

	public:
	static int uTags;
	static bool genTags;
        static int max_step;

	static std::map<int,std::map<int,int>> pTagCnt;
	static std::map<int,std::map<int,std::map<int,std::map<int,std::map<int,int>>>>> tagMap;
	static std::map<int,std::map<int,std::map<int,std::map<int,int>>>> myTagMap;
        static void clearTagMap();
        static void clearMyTagMap();
	static void communicateTags();
        static void registerId(int tid);
	static int tid();
	static volatile int numTeamsFinished;
        static volatile int updateMetadata_request;
        static volatile int updateMetadata_noticed;
        static volatile int updateMetadata_done;
	static Barrier * globalBarrier;
        static void syncProcesses();
	static void multifabBuildFabCon(RegionGraph* graph, const MultiFab& mf, const Periodicity& period);
	static void serviceLocalRequests(RegionGraph *graph, int tg);
	static void serviceRemoteRequests(RegionGraph *graph, int graphID, int nGraphs);
	static void serviceRemoteRequests(RegionGraph *graph);  
	static void serviceSingleGraphComm(RegionGraph* graph, int tid);
	static void serviceMultipleGraphComm(RegionGraph graphArray[], int nGraphs, bool cpyAcross, int tid);
	static void serviceMultipleGraphCommDynamic(std::vector<RegionGraph*> graphArray, bool cpyAcross, int tid);
        static void flattenGraphHierarchy(std::vector<std::vector<RegionGraph*> >graphArray, std::vector<RegionGraph*> &flattenedGraphArray);
	static void serviceMultipleGraphComm(RegionGraph graphArray[], int nGraphs, int tid);
	static void fillBoundaryPush(RegionGraph* graph, MultiFab* mf, int f);
	static void fillBoundaryPull(RegionGraph* graph, MultiFab* mf, int f,  bool singleT);

        void   multifabExtractCopyAssoc(void* threadInfo);
	static void multifabExtractCopyAssoc(RegionGraph* gDst, RegionGraph* gSrc, const MultiFab& dmf, const MultiFab& smf, int nc, int ng, int ngSrc, const Periodicity& period);
	static void multifabExtractCopyAssoc(RegionGraph* gDst, RegionGraph* gSrc, const MultiFab& dmf, const MultiFab& smf, const Periodicity& period);
	static void multifabCopyPushAsync(RegionGraph* destGraph, RegionGraph* srcGraph, MultiFab* dmf, MultiFab* smf, int f, int dstcomp, int srccomp, int nc, int ng, int ngsrc, bool singleT);
	static void multifabCopyPushAsync(RegionGraph* destGraph, RegionGraph* srcGraph, MultiFab* dmf, MultiFab* smf, int f, bool singleT); 
	static void multifabCopyPull(RegionGraph* destGraph, RegionGraph* srcGraph, MultiFab* dmf, MultiFab* smf, int f, int dstcomp, int srccomp, int nc, int ng, int ngsrc, bool singleT);
	static void multifabCopyPull(RegionGraph* destGraph, RegionGraph* srcGraph, MultiFab* dmf, MultiFab* smf, int f, bool singleT);
	static void serviceLocalGridCopyRequests(std::vector<RegionGraph*> graphArray, int g, int tg);
	static void serviceRemoteGridCopyRequests(std::vector<RegionGraph*> graphArray, int g, int nGraph, int tg);
	static void resetRemoteGridCopyRequests(std::vector<RegionGraph*> graphArray, int g, int nGraph, int tg);



        static void multifabCopyPush(RegionGraph* destGraph, RegionGraph* srcGraph, amrex::MultiFab* dmf, amrex::MultiFab* smf, int f, int dstcomp, int srccomp, int nc, int ng, int ngsrc, bool singleT);
        static void multifabCopyPush(RegionGraph* destGraph, RegionGraph* srcGraph, amrex::MultiFab* dmf, amrex::MultiFab* smf, int f, bool singleT);
       
       static void multifabCopyPush_1Team(RegionGraph* destGraph, RegionGraph* srcGraph, amrex::MultiFab* dmf, amrex::MultiFab* smf, int f, int dstcomp, int srccomp, int nc, int ng, int ngsrc, bool singleT);
       static void fillBoundaryPull_1Team(RegionGraph *graph, amrex::MultiFab& mf, int f);

        static void fillBoundaryPush(amrex::RGIter& rgi, amrex::MultiFab& mf);
        static void fillBoundaryPull(amrex::RGIter& rgi, amrex::MultiFab& mf, bool singleT);
        static void fillBoundaryPush(amrex::RGIter& rgi, RegionGraph *graph, amrex::MultiFab& mf);
        static void fillBoundaryPull(amrex::RGIter& rgi, RegionGraph *graph, amrex::MultiFab& mf, bool singleT);

    }; // class Perilla
}
#endif
