
#include <algorithm>

#include <AmrCore.H>
#include <ParmParse.H>

#ifdef USE_PARTICLES
#include <AmrParGDB.H>
#endif

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
{
    Initialize();
    Geometry::Setup();
    int max_level_in = -1;
    Array<int> n_cell_in(BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; i++) n_cell_in[i] = -1;
    InitAmrCore(max_level_in,n_cell_in);    
}

AmrCore::AmrCore (const RealBox* rb, int max_level_in, Array<int> n_cell_in, int coord)
{
    Initialize();
    Geometry::Setup(rb,coord);
    InitAmrCore(max_level_in,n_cell_in);
}

AmrCore::~AmrCore ()
{
#ifdef USE_PARTICLES
    delete m_gdb;
#endif
    Finalize();
}

void
AmrCore::InitAmrCore (int max_level_in, Array<int> n_cell_in)
{
    verbose   = 0;
    grid_eff  = 0.7;
    n_proper  = 1;
    
    ParmParse pp("amr");

    pp.query("v",verbose);

    if (max_level_in == -1) {
       pp.get("max_level", max_level);
    } else {
       max_level = max_level_in;
    }

    int nlev = max_level + 1;
    
    blocking_factor.resize(nlev);
    max_grid_size.resize(nlev);
    n_error_buf.resize(nlev);

    geom.resize(nlev);
    dmap.resize(nlev);
    grids.resize(nlev);
#ifdef USE_PARTICLES
    particle_dmap.resize(nlev);
    particle_grids.resize(nlev);
#endif

    for (int i = 0; i < nlev; ++i) {
	n_error_buf[i] = 1;
        blocking_factor[i] = 8;
        max_grid_size[i] = (BL_SPACEDIM == 2) ? 128 : 32;
    }

    // Make the default ref_ratio = 2 for all levels.
    ref_ratio.resize(max_level);
    for (int i = 0; i < max_level; ++i) {
        ref_ratio[i] = 2 * IntVect::TheUnitVector();
    }

    pp.query("n_proper",n_proper);
    pp.query("grid_eff",grid_eff);
    pp.queryarr("n_error_buf",n_error_buf,0,max_level);

    // Read in the refinement ratio IntVects as integer BL_SPACEDIM-tuples.
    if (max_level > 0)
    {
        const int nratios_vect = max_level*BL_SPACEDIM;

        Array<int> ratios_vect(nratios_vect);

        int got_vect = pp.queryarr("ref_ratio_vect",ratios_vect,0,nratios_vect);

        Array<int> ratios(max_level);

        const int got_int = pp.queryarr("ref_ratio",ratios,0,max_level);
   
        if (got_int == 1 && got_vect == 1 && ParallelDescriptor::IOProcessor())
        {
            BoxLib::Warning("Only input *either* ref_ratio or ref_ratio_vect");
        }
        else if (got_vect == 1)
        {
            int k = 0;
            for (int i = 0; i < max_level; i++)
            {
                for (int n = 0; n < BL_SPACEDIM; n++,k++)
                    ref_ratio[i][n] = ratios_vect[k];
            }
        }
        else if (got_int == 1)
        {
            for (int i = 0; i < max_level; i++)
            {
                for (int n = 0; n < BL_SPACEDIM; n++)
                    ref_ratio[i][n] = ratios[i];
            }
        }
        else
        {
            if (ParallelDescriptor::IOProcessor())
                BoxLib::Warning("Using default ref_ratio = 2 at all levels");
        }
    }

    // Read in max_grid_size.  Use defaults if not explicitly defined.
    int cnt = pp.countval("max_grid_size");
    if (cnt == 1)
    {
        // Set all values to the single available value.
        int the_max_grid_size = 0;
        pp.get("max_grid_size",the_max_grid_size);
        for (int i = 0; i <= max_level; ++i) {
            max_grid_size[i] = the_max_grid_size;
        }
    }
    else if (cnt > 1)
    {
        // Otherwise we expect a vector of max_grid_size values.
        pp.getarr("max_grid_size",max_grid_size,0,max_level+1);
    }

    // Read in the blocking_factors.  Use defaults if not explicitly defined.
    cnt = pp.countval("blocking_factor");
    if (cnt == 1)
    {
        // Set all values to the single available value.
        int the_blocking_factor = 0;
        pp.get("blocking_factor",the_blocking_factor);
        for (int i = 0; i <= max_level; ++i) {
            blocking_factor[i] = the_blocking_factor;
        }
    }
    else if (cnt > 1)
    {
        // Otherwise we expect a vector of blocking factors.
        pp.getarr("blocking_factor",blocking_factor,0,max_level+1);
    }

    // Read computational domain and set geometry.
    {
	Array<int> n_cell(BL_SPACEDIM);
	if (n_cell_in[0] == -1)
	{
	    pp.getarr("n_cell",n_cell,0,BL_SPACEDIM);
	}
	else
	{
	    for (int i = 0; i < BL_SPACEDIM; i++) n_cell[i] = n_cell_in[i];
	}

	IntVect lo(IntVect::TheZeroVector()), hi(n_cell);
	hi -= IntVect::TheUnitVector();
	Box index_domain(lo,hi);
	for (int i = 0; i <= max_level; i++)
	{
	    geom[i].define(index_domain);
	    if (i < max_level)
		index_domain.refine(ref_ratio[i]);
	}

	Real offset[BL_SPACEDIM];
	for (int i = 0; i < BL_SPACEDIM; i++)
	{
	    const Real delta = Geometry::ProbLength(i)/(Real)n_cell[i];
	    offset[i]        = Geometry::ProbLo(i) + delta*lo[i];
	}
	CoordSys::SetOffset(offset);
    }
    
#ifdef USE_PARTICLES
    m_gdb = new AmrParGDB(this);
#endif
}

int
AmrCore::MaxRefRatio (int lev) const
{
    int maxval = 0;
    for (int n = 0; n<BL_SPACEDIM; n++) 
        maxval = std::max(maxval,ref_ratio[lev][n]);
    return maxval;
}

void
AmrCore::MakeDistributionMap (int lev)
{
    dmap[lev] = DistributionMapping(grids[lev], ParallelDescriptor::NProcs());
}


bool
AmrCore::LevelDefined (int lev)
{
    return lev <= max_level && !grids[lev].empty() && !dmap[lev].empty();
}
