
#include <algorithm>

#include <AMReX_AmrCore.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#ifdef USE_PARTICLES
#include <AMReX_AmrParGDB.H>
#endif

namespace amrex {

namespace
{
    bool initialized = false;
}

void
AmrCore::Initialize ()
{
    if (initialized) return;
    initialized = true;
}

void
AmrCore::Finalize ()
{
    initialized = false;
}

AmrCore::AmrCore ()
    : AmrMesh()
{
    Initialize();
    InitAmrCore();
}

AmrCore::AmrCore (const RealBox* rb, int max_level_in, const Array<int>& n_cell_in, int coord)
    : AmrMesh(rb, max_level_in, n_cell_in, coord)
{
    Initialize();
    InitAmrCore();
}

AmrCore::~AmrCore ()
{
    Finalize();
}

void
AmrCore::InitAmrCore ()
{
    verbose   = 0;
    ParmParse pp("amr");
    pp.query("v",verbose);
    
#ifdef USE_PARTICLES
    m_gdb.reset(new AmrParGDB(this));
#endif
}

void
AmrCore::InitFromScratch (Real time)
{
    // define coarse level BoxArray and DistributionMap
    {
	finest_level = 0;

	const BoxArray& ba = MakeBaseGrids();
	DistributionMapping dm(ba);

	MakeNewLevelFromScratch(0, time, ba, dm);

	SetBoxArray(0, ba);
	SetDistributionMap(0, dm);
    }

    if (max_level > 0) // build fine levels
    {
	Array<BoxArray> new_grids(max_level+1);
	new_grids[0] = grids[0];
	do
	{
	    int new_finest;

	    // Add (at most) one level at a time.
	    MakeNewGrids(finest_level,time,new_finest,new_grids);

	    if (new_finest <= finest_level) break;
	    finest_level = new_finest;

	    DistributionMapping dm(new_grids[new_finest]);

	    MakeNewLevelFromScratch(new_finest, time, new_grids[finest_level], dm);

	    SetBoxArray(new_finest, new_grids[new_finest]);
	    SetDistributionMap(new_finest, dm);
	}
	while (finest_level < max_level);

	// Iterate grids to ensure fine grids encompass all interesting junk.
	for (int it=0; it<4; ++it)  // try at most 4 times
	{
	    for (int i = 1; i <= finest_level; ++i) {
		new_grids[i] = grids[i];
	    }

	    int new_finest;
	    MakeNewGrids(0, time, new_finest, new_grids);
	    
	    if (new_finest < finest_level) break;
	    finest_level = new_finest;

	    bool grids_the_same = true;
	    for (int lev = 1; lev <= new_finest; ++lev) {
		if (new_grids[lev] != grids[lev]) {
		    grids_the_same = false;
		    DistributionMapping dm(new_grids[lev]);

		    MakeNewLevelFromScratch(lev, time, new_grids[lev], dm);

		    SetBoxArray(lev, new_grids[lev]);
		    SetDistributionMap(lev, dm);
		}
	    }
	    if (grids_the_same) break;
	}
    }
}

void
AmrCore::regrid (int lbase, Real time, bool)
{
    int new_finest;
    Array<BoxArray> new_grids(finest_level+2);
    MakeNewGrids(lbase, time, new_finest, new_grids);

    BL_ASSERT(new_finest <= finest_level+1);

    for (int lev = lbase+1; lev <= new_finest; ++lev)
    {
	if (lev <= finest_level) // an old level
	{
	    if (new_grids[lev] != grids[lev]) // otherwise nothing
	    {
		DistributionMapping new_dmap(new_grids[lev]);
		RemakeLevel(lev, time, new_grids[lev], new_dmap);
		SetBoxArray(lev, new_grids[lev]);
		SetDistributionMap(lev, new_dmap);
	    }
	}
	else  // a new level
	{
	    DistributionMapping new_dmap(new_grids[lev]);
	    MakeNewLevelFromCoarse(lev, time, new_grids[lev], new_dmap);
	    SetBoxArray(lev, new_grids[lev]);
	    SetDistributionMap(lev, new_dmap);
	}
    }

    for (int lev = new_finest+1; lev <= finest_level; ++lev) {
	ClearLevel(lev);
	ClearBoxArray(lev);
	ClearDistributionMap(lev);
    }

    finest_level = new_finest;
}

}
