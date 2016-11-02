
#include <algorithm>

#include <AmrCore.H>
#include <ParmParse.H>
#include <TagBox.H>
#include <Cluster.H>

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

AmrCore::AmrCore (const RealBox* rb, int max_level_in, const Array<int>& n_cell_in, int coord)
{
    Initialize();
    Geometry::Setup(rb,coord);
    InitAmrCore(max_level_in,n_cell_in);
}

AmrCore::~AmrCore ()
{
    Finalize();
}

void
AmrCore::InitAmrCore (int max_level_in, const Array<int>& n_cell_in)
{
    verbose   = 0;
    grid_eff  = 0.7;
    n_proper  = 1;

    use_fixed_coarse_grids = false;
    use_fixed_upto_level   = 0;

    refine_grid_layout = true;
    
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

    {
	// Are we going to keep the grids at certain levels fixed?
	pp.query("use_fixed_coarse_grids", use_fixed_coarse_grids);

	// chop up grids to have more grids than the number of procs
	pp.query("refine_grid_layout", refine_grid_layout);
    }

    finest_level = -1;
    
#ifdef USE_PARTICLES
    m_gdb = std::unique_ptr<AmrParGDB>(new AmrParGDB(this));
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
AmrCore::SetDistributionMap (int lev, const DistributionMapping& dmap_in)
{ 
    if (dmap[lev] != dmap_in) dmap[lev] = dmap_in;
}

void
AmrCore::SetBoxArray (int lev, const BoxArray& ba_in)
{
    if (grids[lev] != ba_in) grids[lev] = ba_in;
}

void
AmrCore::ClearDistributionMap (int lev)
{ 
    dmap[lev] = DistributionMapping();
}

void
AmrCore::ClearBoxArray (int lev)
{
    grids[lev] = BoxArray();
}

bool
AmrCore::LevelDefined (int lev)
{
    return lev <= max_level && !grids[lev].empty() && !dmap[lev].empty();
}

void
AmrCore::ChopGrids (int lev, BoxArray& ba, int target_size) const
{
    for (int cnt = 1; cnt <= 4; cnt *= 2)
    {
	const int ChunkSize = max_grid_size[lev]/cnt;
	IntVect chunk(D_DECL(ChunkSize,ChunkSize,ChunkSize));

	for (int j = BL_SPACEDIM-1; j >= 0 ; j--)
	{
	    chunk[j] /= 2;
	    
	    if ( (ba.size() < target_size) && (chunk[j]%blocking_factor[lev] == 0) )
	    {
		ba.maxSize(chunk);
	    }
	}
    }
}

BoxArray
AmrCore::MakeBaseGrids () const
{
    BoxArray ba(BoxLib::coarsen(geom[0].Domain(),2));
    ba.maxSize(max_grid_size[0]/2);
    ba.refine(2);
    if (refine_grid_layout) {
	ChopGrids(0, ba, ParallelDescriptor::NProcs());
    }
    if (ba == grids[0]) {
	ba = grids[0];  // to avoid dupliates
    }
    return ba;
}


void
AmrCore::MakeNewGrids (int lbase, Real time, int& new_finest, Array<BoxArray>& new_grids)
{
    BL_ASSERT(lbase < max_level);

    // Add at most one new level
    int max_crse = std::min(finest_level, max_level-1);

    if (new_grids.size() < max_crse+2) new_grids.resize(max_crse+2);

    //
    // Construct problem domain at each level.
    //
    Array<IntVect> bf_lev(max_level); // Blocking factor at each level.
    Array<IntVect> rr_lev(max_level);
    Array<Box>     pc_domain(max_level);  // Coarsened problem domain.

    for (int i = 0; i <= max_crse; i++)
    {
        for (int n=0; n<BL_SPACEDIM; n++)
            bf_lev[i][n] = std::max(1,blocking_factor[i+1]/ref_ratio[i][n]);
    }
    for (int i = lbase; i < max_crse; i++)
    {
        for (int n=0; n<BL_SPACEDIM; n++)
            rr_lev[i][n] = (ref_ratio[i][n]*bf_lev[i][n])/bf_lev[i+1][n];
    }
    for (int i = lbase; i <= max_crse; i++) {
	pc_domain[i] = BoxLib::coarsen(Geom(i).Domain(),bf_lev[i]);
    }
    //
    // Construct proper nesting domains.
    //
    Array<BoxList> p_n(max_level);      // Proper nesting domain.
    Array<BoxList> p_n_comp(max_level); // Complement proper nesting domain.

    BoxList bl(grids[lbase]);
    bl.simplify();
    bl.coarsen(bf_lev[lbase]);
    p_n_comp[lbase].complementIn(pc_domain[lbase],bl);
    p_n_comp[lbase].simplify();
    p_n_comp[lbase].accrete(n_proper);
    if (Geometry::isAnyPeriodic()) {
	ProjPeriodic(p_n_comp[lbase], Geometry(pc_domain[lbase]));
    }
    p_n[lbase].complementIn(pc_domain[lbase],p_n_comp[lbase]);
    p_n[lbase].simplify();
    bl.clear();

    for (int i = lbase+1; i <= max_crse; i++)
    {
        p_n_comp[i] = p_n_comp[i-1];

        // Need to simplify p_n_comp or the number of grids can too large for many levels.
        p_n_comp[i].simplify();

        p_n_comp[i].refine(rr_lev[i-1]);
        p_n_comp[i].accrete(n_proper);

	if (Geometry::isAnyPeriodic()) {
	    ProjPeriodic(p_n_comp[i], Geometry(pc_domain[i]));
	}

        p_n[i].complementIn(pc_domain[i],p_n_comp[i]);
        p_n[i].simplify();
    }
    //
    // Now generate grids from finest level down.
    //
    new_finest = lbase;

    for (int levc = max_crse; levc >= lbase; levc--)
    {
        int levf = levc+1;
        //
        // Construct TagBoxArray with sufficient grow factor to contain
        // new levels projected down to this level.
        //
        int ngrow = 0;

        if (levf < new_finest)
        {
            BoxArray ba_proj(new_grids[levf+1]);

            ba_proj.coarsen(ref_ratio[levf]);
            ba_proj.growcoarsen(n_proper, ref_ratio[levc]);

            BoxArray levcBA = grids[levc];

            while (!levcBA.contains(ba_proj))
            {
                levcBA.grow(1);
                ++ngrow;
            }
        }
        TagBoxArray tags(grids[levc],n_error_buf[levc]+ngrow);
    
        //
        // Only use error estimation to tag cells for the creation of new grids
        //      if the grids at that level aren't already fixed.
        //

        if ( ! (useFixedCoarseGrids() && levc < useFixedUpToLevel()) ) {
	    ErrorEst(levc, tags, time, ngrow);
	}

        //
        // If new grids have been constructed above this level, project
        // those grids down and tag cells on intersections to ensure
        // proper nesting.
        //
        // NOTE: this loop replaces the previous code:
        //      if (levf < new_finest) 
        //          tags.setVal(ba_proj,TagBox::SET);
        // The problem with this code is that it effectively 
        // "buffered the buffer cells",  i.e., the grids at level
        // levf+1 which were created by buffering with n_error_buf[levf]
        // are then coarsened down twice to define tagging at
        // level levc, which will then also be buffered.  This can
        // create grids which are larger than necessary.
        //
        if (levf < new_finest)
        {
            int nerr = n_error_buf[levf];

            BoxList bl_tagged(new_grids[levf+1]);
            bl_tagged.simplify();
            bl_tagged.coarsen(ref_ratio[levf]);
            //
            // This grows the boxes by nerr if they touch the edge of the
            // domain in preparation for them being shrunk by nerr later.
            // We want the net effect to be that grids are NOT shrunk away
            // from the edges of the domain.
            //
            for (BoxList::iterator blt = bl_tagged.begin(), End = bl_tagged.end();
                 blt != End;
                 ++blt)
            {
                for (int idir = 0; idir < BL_SPACEDIM; idir++)
                {
                    if (blt->smallEnd(idir) == Geom(levf).Domain().smallEnd(idir))
                        blt->growLo(idir,nerr);
                    if (blt->bigEnd(idir) == Geom(levf).Domain().bigEnd(idir))
                        blt->growHi(idir,nerr);
                }
            }
            Box mboxF = BoxLib::grow(bl_tagged.minimalBox(),1);
            BoxList blFcomp;
            blFcomp.complementIn(mboxF,bl_tagged);
            blFcomp.simplify();
            bl_tagged.clear();

            const IntVect& iv = IntVect(D_DECL(nerr/ref_ratio[levf][0],
                                               nerr/ref_ratio[levf][1],
                                               nerr/ref_ratio[levf][2]));
            blFcomp.accrete(iv);
            BoxList blF;
            blF.complementIn(mboxF,blFcomp);
            BoxArray baF(blF);
            blF.clear();
            baF.grow(n_proper);
            //
            // We need to do this in case the error buffering at
            // levc will not be enough to cover the error buffering
            // at levf which was just subtracted off.
            //
            for (int idir = 0; idir < BL_SPACEDIM; idir++) 
            {
                if (nerr > n_error_buf[levc]*ref_ratio[levc][idir]) 
                    baF.grow(idir,nerr-n_error_buf[levc]*ref_ratio[levc][idir]);
            }

            baF.coarsen(ref_ratio[levc]);

            tags.setVal(baF,TagBox::SET);
        }
        //
        // Buffer error cells.
        //
        tags.buffer(n_error_buf[levc]+ngrow);

        if (useFixedCoarseGrids())
        {
	    if (levc>=useFixedUpToLevel())
	    {
		tags.setVal(GetAreaNotToTag(levc), TagBox::CLEAR);
	    }
	    else
	    {
		new_finest = std::max(new_finest,levf);
	    }
        }

        //
        // Coarsen the taglist by blocking_factor/ref_ratio.
        //
        int bl_max = 0;
        for (int n=0; n<BL_SPACEDIM; n++)
            bl_max = std::max(bl_max,bf_lev[levc][n]);
        if (bl_max > 1) 
            tags.coarsen(bf_lev[levc]);
        //
        // Remove or add tagged points which violate/satisfy additional 
        // user-specified criteria.
        //
	ManualTagsPlacement(levc, tags, bf_lev);
        //
        // Map tagged points through periodic boundaries, if any.
        //
        tags.mapPeriodic(Geometry(pc_domain[levc]));
        //
        // Remove cells outside proper nesting domain for this level.
        //
        tags.setVal(p_n_comp[levc],TagBox::CLEAR);
        //
        // Create initial cluster containing all tagged points.
        //
	std::vector<IntVect> tagvec;
	tags.collate(tagvec);
        tags.clear();

        if (tagvec.size() > 0)
        {
            //
            // Created new level, now generate efficient grids.
            //
            if ( !(useFixedCoarseGrids() && levc<useFixedUpToLevel()) ) {
                new_finest = std::max(new_finest,levf);
	    }
            //
            // Construct initial cluster.
            //
            ClusterList clist(&tagvec[0], tagvec.size());
            clist.chop(grid_eff);
            BoxDomain bd;
            bd.add(p_n[levc]);
            clist.intersect(bd);
            bd.clear();
            //
            // Efficient properly nested Clusters have been constructed
            // now generate list of grids at level levf.
            //
            BoxList new_bx;
            clist.boxList(new_bx);
            new_bx.refine(bf_lev[levc]);
            new_bx.simplify();
            BL_ASSERT(new_bx.isDisjoint());
	    
	    if (new_bx.size()>0) {
		if ( !(Geom(levc).Domain().contains(BoxArray(new_bx).minimalBox())) ) {
		// Chop new grids outside domain, note that this is likely to result in
		//  new grids that violate blocking_factor....see warning checking below
		    new_bx = BoxLib::intersect(new_bx,Geom(levc).Domain());
		}
	    }

            IntVect largest_grid_size;
            for (int n = 0; n < BL_SPACEDIM; n++)
                largest_grid_size[n] = max_grid_size[levf] / ref_ratio[levc][n];
            //
            // Ensure new grid boxes are at most max_grid_size in index dirs.
            //
            new_bx.maxSize(largest_grid_size);

#ifdef BL_FIX_GATHERV_ERROR
	      int wcount = 0, iLGS = largest_grid_size[0];

              while (new_bx.size() < 64 && wcount++ < 4)
              {
                  iLGS /= 2;
                  if (ParallelDescriptor::IOProcessor())
                  {
                      std::cout << "BL_FIX_GATHERV_ERROR:  using iLGS = " << iLGS
                                << "   largest_grid_size was:  " << largest_grid_size[0]
                                << '\n';
                      std::cout << "BL_FIX_GATHERV_ERROR:  new_bx.size() was:   "
                                << new_bx.size() << '\n';
                  }

                  new_bx.maxSize(iLGS);

                  if (ParallelDescriptor::IOProcessor())
                  {
                      std::cout << "BL_FIX_GATHERV_ERROR:  new_bx.size() now:   "
                                << new_bx.size() << '\n';
                  }
	      }
#endif
            //
            // Refine up to levf.
            //
            new_bx.refine(ref_ratio[levc]);
            BL_ASSERT(new_bx.isDisjoint());

	    if (new_bx.size()>0) {
		if ( !(Geom(levf).Domain().contains(BoxArray(new_bx).minimalBox())) ) {
		    new_bx = BoxLib::intersect(new_bx,Geom(levf).Domain());
		}
		if (ParallelDescriptor::IOProcessor()) {
		    for (int d=0; d<BL_SPACEDIM; ++d) {
			bool ok = true;
			for (BoxList::const_iterator bli = new_bx.begin(); bli != new_bx.end(); ++bli) {
			    int len = bli->length(d);
			    int bf = blocking_factor[levf];
			    ok &= (len/bf) * bf == len;
			}
			if (!ok) {
			    BoxLib::Warning("WARNING: New grids violate blocking factor near upper boundary");
			}
		    }
		}
	    }

            if(levf > useFixedUpToLevel()) {
              new_grids[levf].define(new_bx);
	    }
        }
    }

    if (refine_grid_layout) {
	for (int lev = lbase+1; lev <= new_finest; ++lev) {
	    ChopGrids(lev,new_grids[lev],ParallelDescriptor::NProcs());
	    if (new_grids[lev] == grids[lev]) {
		new_grids[lev] = grids[lev]; // to avoid dupliates
	    }
	}
    }
}

void
AmrCore::ProjPeriodic (BoxList& blout, const Geometry& geom)
{
    //
    // Add periodic translates to blout.
    //
    Box domain = geom.Domain();

    BoxList blorig(blout);

    int nist,njst,nkst;
    int niend,njend,nkend;
    nist = njst = nkst = 0;
    niend = njend = nkend = 0;
    D_TERM( nist , =njst , =nkst ) = -1;
    D_TERM( niend , =njend , =nkend ) = +1;

    int ri,rj,rk;
    for (ri = nist; ri <= niend; ri++)
    {
        if (ri != 0 && !geom.isPeriodic(0))
            continue;
        if (ri != 0 && geom.isPeriodic(0))
            blorig.shift(0,ri*domain.length(0));
        for (rj = njst; rj <= njend; rj++)
        {
            if (rj != 0 && !geom.isPeriodic(1))
                continue;
            if (rj != 0 && geom.isPeriodic(1))
                blorig.shift(1,rj*domain.length(1));
            for (rk = nkst; rk <= nkend; rk++)
            {
                if (rk != 0 && !geom.isPeriodic(2))
                    continue;
                if (rk != 0 && geom.isPeriodic(2))
                    blorig.shift(2,rk*domain.length(2));

                BoxList tmp(blorig);
                tmp.intersect(domain);
                blout.catenate(tmp);
 
                if (rk != 0 && geom.isPeriodic(2))
                    blorig.shift(2,-rk*domain.length(2));
            }
            if (rj != 0 && geom.isPeriodic(1))
                blorig.shift(1,-rj*domain.length(1));
        }
        if (ri != 0 && geom.isPeriodic(0))
            blorig.shift(0,-ri*domain.length(0));
    }
}
