
#include <algorithm>

#include <AMReX_AmrCore.H>
#include <AMReX_Print.H>

#ifdef AMREX_PARTICLES
#include <AMReX_AmrParGDB.H>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

AmrCore::AmrCore ()
    : AmrMesh()
{
    InitAmrCore();
}

AmrCore::AmrCore (const RealBox* rb, int max_level_in,
                  const Vector<int>& n_cell_in, int coord,
                  Vector<IntVect> ref_ratios, const int* is_per)
    : AmrMesh(rb, max_level_in, n_cell_in, coord, std::move(ref_ratios), is_per)
{
    InitAmrCore();
}

AmrCore::AmrCore (const RealBox& rb, int max_level_in,
                  const Vector<int>& n_cell_in, int coord,
                  Vector<IntVect> const& ref_ratios,
                  Array<int,AMREX_SPACEDIM> const& is_per)
    : AmrMesh(rb, max_level_in, n_cell_in, coord, ref_ratios, is_per)
{
    InitAmrCore();
}

AmrCore::AmrCore (Geometry const& level_0_gome, AmrInfo const& amr_info)
    : AmrMesh(level_0_gome,amr_info)
{
#ifdef AMREX_PARTICLES
    m_gdb.reset(new AmrParGDB(this));
#endif
}

AmrCore::~AmrCore ()
{
}

void
AmrCore::InitAmrCore ()
{
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
    if (lbase >= max_level) return;

    int new_finest;
    Vector<BoxArray> new_grids(finest_level+2);
    MakeNewGrids(lbase, time, new_finest, new_grids);

    BL_ASSERT(new_finest <= finest_level+1);

    bool coarse_ba_changed = false;
    for (int lev = lbase+1; lev <= new_finest; ++lev)
    {
        if (lev <= finest_level) // an old level
        {
            bool ba_changed = (new_grids[lev] != grids[lev]);
	    if (ba_changed or coarse_ba_changed) {
                BoxArray level_grids = grids[lev];
                DistributionMapping level_dmap = dmap[lev];
                if (ba_changed) {
                    level_grids = new_grids[lev];
                    level_dmap = DistributionMapping(level_grids);
                }
                const auto old_num_setdm = num_setdm;
                RemakeLevel(lev, time, level_grids, level_dmap);
                SetBoxArray(lev, level_grids);
                if (old_num_setdm == num_setdm) {
                    SetDistributionMap(lev, level_dmap);
                }
            }
            coarse_ba_changed = ba_changed;;
	}
	else  // a new level
	{
            DistributionMapping new_dmap(new_grids[lev]);
            const auto old_num_setdm = num_setdm;
            MakeNewLevelFromCoarse(lev, time, new_grids[lev], new_dmap);
            SetBoxArray(lev, new_grids[lev]);
            if (old_num_setdm == num_setdm) {
                SetDistributionMap(lev, new_dmap);
            }
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
        Long                      ncells  = bs.numPts();
        double                    ntot    = Geom(lev).Domain().d_numPts();
        Real                      frac    = Real(100.0)*(Real(ncells) / ntot);

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
	    Long vmin = std::numeric_limits<Long>::max();
	    Long vmax = -1;
	    int lmax = -1;
	    int smin = std::numeric_limits<int>::max();
            int imax = std::numeric_limits<int>::lowest();
            int imin = std::numeric_limits<int>::lowest();
#ifdef _OPENMP
#pragma omp parallel
#endif	    
	    {
		Long vmin_this = std::numeric_limits<Long>::max();
		Long vmax_this = -1;
		int lmax_this = -1;
		int smin_this = std::numeric_limits<int>::max();
                int imax_this = std::numeric_limits<int>::lowest();
                int imin_this = std::numeric_limits<int>::lowest();
#ifdef _OPENMP
#pragma omp for
#endif	    	    
		for (int k = 0; k < numgrid; k++) {
		    const Box& bx = bs[k];
		    Long v = bx.volume();
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
