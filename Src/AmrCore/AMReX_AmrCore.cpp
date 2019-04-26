
#include <algorithm>

#include <AMReX_AmrCore.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#ifdef AMREX_PARTICLES
#include <AMReX_AmrParGDB.H>
#endif

#ifdef _OPENMP
#include <omp.h>
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

AmrCore::AmrCore (const RealBox* rb, int max_level_in, const Vector<int>& n_cell_in, int coord, Vector<IntVect> ref_ratios)
  : AmrMesh(rb, max_level_in, n_cell_in, coord, std::move(ref_ratios))
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

#ifdef AMREX_PARTICLES
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
    Vector<BoxArray> new_grids(finest_level+2);
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


void
AmrCore::printGridSummary (std::ostream& os, int min_lev, int max_lev) const noexcept
{
    for (int lev = min_lev; lev <= max_lev; lev++)
    {
        const BoxArray&           bs      = boxArray(lev);
        int                       numgrid = bs.size();
        long                      ncells  = bs.numPts();
        double                    ntot    = Geom(lev).Domain().d_numPts();
        Real                      frac    = 100.0*(Real(ncells) / ntot);

        os << "  Level "
           << lev
           << "   "
           << numgrid
           << " grids  "
           << ncells
           << " cells  "
           << frac
           << " % of domain"
           << '\n';

	if (numgrid > 1) {
	    long vmin = std::numeric_limits<long>::max();
	    long vmax = -1;
	    int lmax = -1;
	    int smin = std::numeric_limits<int>::max();
	    int imax, imin;
#ifdef _OPENMP
#pragma omp parallel
#endif	    
	    {
		long vmin_this = std::numeric_limits<long>::max();
		long vmax_this = -1;
		int lmax_this = -1;
		int smin_this = std::numeric_limits<int>::max();
		int imax_this, imin_this;
#ifdef _OPENMP
#pragma omp for
#endif	    	    
		for (int k = 0; k < numgrid; k++) {
		    const Box& bx = bs[k];
		    long v = bx.volume();
		    int ss = bx.shortside();
		    int ls = bx.longside();
		    if (v < vmin_this || (v == vmin_this && ss < smin_this)) {
			vmin_this = v;
			smin_this = ss;
			imin_this = k;
		    }
		    if (v > vmax_this || (v == vmax_this && ls > lmax_this)) {
			vmax_this = v;
			lmax_this = ls;
			imax_this = k;
		    }
		}
#ifdef _OPENMP
#pragma omp critical (amr_prtgs)
#endif	    	    
		{
		    if (vmin_this < vmin || (vmin_this == vmin && smin_this < smin)) {
			vmin = vmin_this;
			smin = smin_this;
			imin = imin_this;
		    }
		    if (vmax_this > vmax || (vmax_this == vmax && lmax_this > lmax)) {
			vmax = vmax_this;
			lmax = lmax_this;
			imax = imax_this;
		    }
		}
	    }
	    const Box& bmin = bs[imin];
	    const Box& bmax = bs[imax];
	    os << "           "
	       << " smallest grid: "
		AMREX_D_TERM(<< bmin.length(0),
		       << " x " << bmin.length(1),
		       << " x " << bmin.length(2))
	       << "  biggest grid: "
		AMREX_D_TERM(<< bmax.length(0),
		       << " x " << bmax.length(1),
		       << " x " << bmax.length(2))
	       << '\n';
	}
    }

    os << std::endl; // Make sure we flush!
}

}
