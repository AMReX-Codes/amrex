
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
    MakeNewGrids(time);
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
