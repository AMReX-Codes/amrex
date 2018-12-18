#if 0
#ifndef RGITER_H_
#define RGITER_H_

#include <RegionGraph.H>
#include <fstream>

namespace amrex{

    class AsyncFillPatchIterator;

    class RGIter
    {
	public:
	    int tid;
	    int ntid;
	    int tg;
	    int currentRegion;
	    int currentTile;
	    int totalItr;
	    int currentItr;
	    bool tiling;
	    bool implicit;
	    bool ppteams;
	    bool haveDepGraph;
	    RegionGraph* itrGraph;
	    RegionGraph* depGraph;
	    int  boxGrow, index, scomp, ncomp, iteration;
	    double time;
	    double getFireableTime;

	    IndexType     typ;

	    std::vector<amrex::AsyncFillPatchIterator*> m_level_afpi; 
	    std::vector<amrex::AsyncFillPatchIterator*> m_upper_level_afpi;
	    std::ofstream fout;

	public:
	    RGIter(RegionGraph* rg, bool enableAllTasks=false);
	    RGIter(RegionGraph* rg, RegionGraph* drg, bool isDep=true);
	    RGIter(amrex::AsyncFillPatchIterator* afpi, bool enableAllTasks=false);
	    RGIter(std::vector<amrex::AsyncFillPatchIterator*> afpi, std::vector<amrex::AsyncFillPatchIterator*> upper_afpi, 
		    amrex::MultiFab& dest, int  bG, double tm, int  ind, int  sc, int nc, int itr);
	    ~RGIter();

	    void init();
	    void sync_workers();
	    //! Increment iterator to the next tile we own.
	    void operator++ ();
	    //! Is the iterator valid, are more regions to iterate over?
	    bool isValid();
	    int LocalIndex() const { return currentRegion; }

	    amrex::Box tileBox();
	    amrex::Box validBox() const;
	    amrex::Box tilebox();
	    amrex::Box growntilebox();
	    amrex::Box growntilebox(int ng);
	    amrex::Box nodaltilebox(int dir);
    };
}
#endif
#endif
