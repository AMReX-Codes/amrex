
#include <AMReX.H>
#include <AMReX_AmrMesh.H>
#include <AMReX_Cluster.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

namespace amrex {

namespace
{
    bool initialized = false;
}

void
AmrMesh::Initialize ()
{
    if (initialized) return;
    initialized = true;
}

void
AmrMesh::Finalize ()
{
    initialized = false;
}


AmrMesh::AmrMesh ()
{
    Initialize();
    Geometry::Setup();
    int max_level_in = -1;
    Vector<int> n_cell_in(AMREX_SPACEDIM, -1);
    InitAmrMesh(max_level_in,n_cell_in);
}

AmrMesh::AmrMesh (const RealBox* rb, int max_level_in, const Vector<int>& n_cell_in, int coord,
                  Vector<IntVect> a_refrat)
{
  Initialize();

  Geometry::Setup(rb,coord);
  InitAmrMesh(max_level_in,n_cell_in, std::move(a_refrat));
}

AmrMesh::~AmrMesh ()
{
    Finalize();
}

void
AmrMesh::InitAmrMesh (int max_level_in, const Vector<int>& n_cell_in, Vector<IntVect> a_refrat)
{
    verbose   = 0;
    grid_eff  = 0.7;
    n_proper  = 1;

    use_fixed_coarse_grids = false;
    use_fixed_upto_level   = 0;
    refine_grid_layout     = true;
    check_input            = true;

    use_new_chop         = false;
    iterate_on_new_grids = true;

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
        n_error_buf[i]     = IntVect{AMREX_D_DECL(1,1,1)};
        blocking_factor[i] = IntVect{AMREX_D_DECL(8,8,8)};
        max_grid_size[i]   = (AMREX_SPACEDIM == 2) ? IntVect{AMREX_D_DECL(128,128,128)}
                                                   : IntVect{AMREX_D_DECL(32,32,32)};
    }

    // Make the default ref_ratio = 2 for all levels.
    ref_ratio.resize(max_level);
    for (int i = 0; i < max_level; ++i)
    {
      ref_ratio[i] = 2 * IntVect::TheUnitVector();
    }

    pp.query("n_proper",n_proper);
    pp.query("grid_eff",grid_eff);
    int cnt = pp.countval("n_error_buf");
    if (cnt > 0) {
        Vector<int> neb;
        pp.getarr("n_error_buf",neb);
        int n = std::min(cnt, max_level+1);
        for (int i = 0; i < n; ++i) {
            n_error_buf[i] = IntVect(neb[i]);
        }
        for (int i = n; i <= max_level; ++i) {
            n_error_buf[i] = IntVect(neb[cnt-1]);
        }
    }

    cnt = pp.countval("n_error_buf_x");
    if (cnt > 0) {
        int idim = 0;
        Vector<int> neb;
        pp.getarr("n_error_buf_x",neb);
        int n = std::min(cnt, max_level+1);
        for (int i = 0; i < n; ++i) {
            n_error_buf[i][idim] = neb[i];
        }
        for (int i = n; i <= max_level; ++i) {
            n_error_buf[i][idim] = neb[n-1];
        }
    }

#if (AMREX_SPACEDIM > 1)
    cnt = pp.countval("n_error_buf_y");
    if (cnt > 0) {
        int idim = 1;
        Vector<int> neb;
        pp.getarr("n_error_buf_y",neb);
        int n = std::min(cnt, max_level+1);
        for (int i = 0; i < n; ++i) {
            n_error_buf[i][idim] = neb[i];
        }
        for (int i = n; i <= max_level; ++i) {
            n_error_buf[i][idim] = neb[n-1];
        }
    }
#endif

#if (AMREX_SPACEDIM == 3)
    cnt = pp.countval("n_error_buf_z");
    if (cnt > 0) {
        int idim = 2;
        Vector<int> neb;
        pp.getarr("n_error_buf_z",neb);
        int n = std::min(cnt, max_level+1);
        for (int i = 0; i < n; ++i) {
            n_error_buf[i][idim] = neb[i];
        }
        for (int i = n; i <= max_level; ++i) {
            n_error_buf[i][idim] = neb[n-1];
        }
    }
#endif

    // Read in the refinement ratio IntVects as integer AMREX_SPACEDIM-tuples.
    if (max_level > 0)
    {
        const int nratios_vect = max_level*AMREX_SPACEDIM;

        Vector<int> ratios_vect(nratios_vect);

        int got_vect = pp.queryarr("ref_ratio_vect",ratios_vect,0,nratios_vect);

        Vector<int> ratios;

        const int got_int = pp.queryarr("ref_ratio",ratios);

        if (got_int == 1 && got_vect == 1)
        {
            amrex::Abort("Only input *either* ref_ratio or ref_ratio_vect");
        }
        else if (got_vect == 1)
        {
            int k = 0;
            for (int i = 0; i < max_level; i++)
            {
                for (int n = 0; n < AMREX_SPACEDIM; n++,k++)
                    ref_ratio[i][n] = ratios_vect[k];
            }
        }
        else if (got_int == 1)
        {
            const int ncnt = ratios.size();
            for (int i = 0; i < ncnt && i < max_level; ++i)
            {
                for (int n = 0; n < AMREX_SPACEDIM; n++) {
                    ref_ratio[i][n] = ratios[i];
                }
            }
            for (int i = ncnt; i < max_level; ++i)
            {
                for (int n = 0; n < AMREX_SPACEDIM; n++) {
                    ref_ratio[i][n] = ratios.back();
                }
            }
        }
        else
        {
            if (verbose) {
                amrex::Print() << "Using default ref_ratio = 2 at all levels\n";
            }
        }
    }
    //if sent in, this wins over everything.
    if(a_refrat.size() > 0)
    {
      for (int i = 0; i < max_level; i++)
      {
          ref_ratio[i] = a_refrat[i];
      }
    }

    // Read in max_grid_size.  Use defaults if not explicitly defined.
    cnt = pp.countval("max_grid_size");
    if (cnt > 0) {
        Vector<int> mgs;
        pp.getarr("max_grid_size",mgs);
        int last_mgs = mgs.back();
        mgs.resize(max_level+1,last_mgs);
        SetMaxGridSize(mgs);
    }

    cnt = pp.countval("max_grid_size_x");
    if (cnt > 0) {
        int idim = 0;
        Vector<int> mgs;
        pp.getarr("max_grid_size_x",mgs);
        int n = std::min(cnt, max_level+1);
        for (int i = 0; i < n; ++i) {
            max_grid_size[i][idim] = mgs[i];
        }
        for (int i = n; i <= max_level; ++i) {
            max_grid_size[i][idim] = mgs[n-1];
        }
    }

#if (AMREX_SPACEDIM > 1)
    cnt = pp.countval("max_grid_size_y");
    if (cnt > 0) {
        int idim = 1;
        Vector<int> mgs;
        pp.getarr("max_grid_size_y",mgs);
        int n = std::min(cnt, max_level+1);
        for (int i = 0; i < n; ++i) {
            max_grid_size[i][idim] = mgs[i];
        }
        for (int i = n; i <= max_level; ++i) {
            max_grid_size[i][idim] = mgs[n-1];
        }
    }
#endif

#if (AMREX_SPACEDIM == 3)
    cnt = pp.countval("max_grid_size_z");
    if (cnt > 0) {
        int idim = 2;
        Vector<int> mgs;
        pp.getarr("max_grid_size_z",mgs);
        int n = std::min(cnt, max_level+1);
        for (int i = 0; i < n; ++i) {
            max_grid_size[i][idim] = mgs[i];
        }
        for (int i = n; i <= max_level; ++i) {
            max_grid_size[i][idim] = mgs[n-1];
        }
    }
#endif

    // Read in the blocking_factors.  Use defaults if not explicitly defined.
    cnt = pp.countval("blocking_factor");
    if (cnt > 0) {
        Vector<int> bf;
        pp.getarr("blocking_factor",bf);
        int last_bf = bf.back();
        bf.resize(max_level+1,last_bf);
        SetBlockingFactor(bf);
    }

    cnt = pp.countval("blocking_factor_x");
    if (cnt > 0) {
        int idim = 0;
        Vector<int> bf;
        pp.getarr("blocking_factor_x",bf);
        int n = std::min(cnt, max_level+1);
        for (int i = 0; i < n; ++i) {
            blocking_factor[i][idim] = bf[i];
        }
        for (int i = n; i <= max_level; ++i) {
            blocking_factor[i][idim] = bf[n-1];
        }
    }

#if (AMREX_SPACEDIM > 1)
    cnt = pp.countval("blocking_factor_y");
    if (cnt > 0) {
        int idim = 1;
        Vector<int> bf;
        pp.getarr("blocking_factor_y",bf);
        int n = std::min(cnt, max_level+1);
        for (int i = 0; i < n; ++i) {
            blocking_factor[i][idim] = bf[i];
        }
        for (int i = n; i <= max_level; ++i) {
            blocking_factor[i][idim] = bf[n-1];
        }
    }
#endif

#if (AMREX_SPACEDIM == 3)
    cnt = pp.countval("blocking_factor_z");
    if (cnt > 0) {
        int idim = 2;
        Vector<int> bf;
        pp.getarr("blocking_factor_z",bf);
        int n = std::min(cnt, max_level+1);
        for (int i = 0; i < n; ++i) {
            blocking_factor[i][idim] = bf[i];
        }
        for (int i = n; i <= max_level; ++i) {
            blocking_factor[i][idim] = bf[n-1];
        }
    }
#endif

    // Read computational domain and set geometry.
    {
	Vector<int> n_cell(AMREX_SPACEDIM);
	if (n_cell_in[0] == -1)
	{
	    pp.getarr("n_cell",n_cell,0,AMREX_SPACEDIM);
	}
	else
	{
	    for (int i = 0; i < AMREX_SPACEDIM; i++) n_cell[i] = n_cell_in[i];
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

	Real offset[AMREX_SPACEDIM];
	for (int i = 0; i < AMREX_SPACEDIM; i++)
	{
	    const Real delta = Geometry::ProbLength(i)/(Real)n_cell[i];
	    offset[i]        = Geometry::ProbLo(i) + delta*lo[i];
	}
	CoordSys::SetOffset(offset);
    }

    {
	// chop up grids to have more grids than the number of procs
	pp.query("refine_grid_layout", refine_grid_layout);
    }

    pp.query("check_input", check_input);

    finest_level = -1;

    if (check_input) checkInput();
}

int
AmrMesh::MaxRefRatio (int lev) const noexcept
{
    int maxval = 0;
    for (int n = 0; n<AMREX_SPACEDIM; n++)
        maxval = std::max(maxval,ref_ratio[lev][n]);
    return maxval;
}

void
AmrMesh::SetDistributionMap (int lev, const DistributionMapping& dmap_in) noexcept
{
    if (dmap[lev] != dmap_in) dmap[lev] = dmap_in;
}

void
AmrMesh::SetBoxArray (int lev, const BoxArray& ba_in) noexcept
{
    if (grids[lev] != ba_in) grids[lev] = ba_in;
}

void
AmrMesh::ClearDistributionMap (int lev) noexcept
{
    dmap[lev] = DistributionMapping();
}

void
AmrMesh::ClearBoxArray (int lev) noexcept
{
    grids[lev] = BoxArray();
}

bool
AmrMesh::LevelDefined (int lev) noexcept
{
    return lev <= max_level && !grids[lev].empty() && !dmap[lev].empty();
}

void
AmrMesh::ChopGrids (int lev, BoxArray& ba, int target_size) const
{
    for (int cnt = 1; cnt <= 4; cnt *= 2)
    {
        IntVect chunk = max_grid_size[lev] / cnt;

	for (int j = AMREX_SPACEDIM-1; j >= 0 ; j--)
	{
	    chunk[j] /= 2;

	    if ( (ba.size() < target_size) && (chunk[j]%blocking_factor[lev][j] == 0) )
	    {
		ba.maxSize(chunk);
	    }
	}
    }
}

BoxArray
AmrMesh::MakeBaseGrids () const
{
    IntVect fac(2);
    const Box& dom = geom[0].Domain();
    const Box dom2 = amrex::refine(amrex::coarsen(dom,2),2);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (dom.length(idim) != dom2.length(idim)) {
            fac[idim] = 1;
        }
    }
    BoxArray ba(amrex::coarsen(dom,fac));
    ba.maxSize(max_grid_size[0]/fac);
    ba.refine(fac);
    // Boxes in ba have even number of cells in each direction
    // unless the domain has odd number of cells in that direction.
    if (refine_grid_layout) {
        ChopGrids(0, ba, ParallelDescriptor::NProcs());
    }
    if (ba == grids[0]) {
        ba = grids[0];  // to avoid duplicates
    }
    return ba;
}


void
AmrMesh::MakeNewGrids (int lbase, Real time, int& new_finest, Vector<BoxArray>& new_grids)
{
    BL_PROFILE("AmrMesh::MakeNewGrids()");

    BL_ASSERT(lbase < max_level);

    // Add at most one new level
    int max_crse = std::min(finest_level, max_level-1);

    if (new_grids.size() < max_crse+2) new_grids.resize(max_crse+2);

    //
    // Construct problem domain at each level.
    //
    Vector<IntVect> bf_lev(max_level); // Blocking factor at each level.
    Vector<IntVect> rr_lev(max_level);
    Vector<Box>     pc_domain(max_level);  // Coarsened problem domain.

    for (int i = 0; i <= max_crse; i++)
    {
        for (int n=0; n<AMREX_SPACEDIM; n++) {
            bf_lev[i][n] = std::max(1,blocking_factor[i+1][n]/ref_ratio[i][n]);
        }
    }
    for (int i = lbase; i < max_crse; i++)
    {
        for (int n=0; n<AMREX_SPACEDIM; n++)
            rr_lev[i][n] = (ref_ratio[i][n]*bf_lev[i][n])/bf_lev[i+1][n];
    }
    for (int i = lbase; i <= max_crse; i++) {
	pc_domain[i] = amrex::coarsen(Geom(i).Domain(),bf_lev[i]);
    }
    //
    // Construct proper nesting domains.
    //
    Vector<BoxList> p_n(max_level);      // Proper nesting domain.
    Vector<BoxList> p_n_comp(max_level); // Complement proper nesting domain.

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
        TagBoxArray tags(grids[levc],dmap[levc],n_error_buf[levc]+ngrow);

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
        // levf+1 which were created by buffering with n_error_buf[levf][idim]
        // are then coarsened down twice to define tagging at
        // level levc, which will then also be buffered.  This can
        // create grids which are larger than necessary.
        //
        if (levf < new_finest)
        {
            // Replace this by n_error_buf that may be anisotropic
            // int nerr = n_error_buf[levf];

            BoxList bl_tagged(new_grids[levf+1]);
            bl_tagged.simplify();
            bl_tagged.coarsen(ref_ratio[levf]);
            //
            // This grows the boxes by n_error_buf[levf][idir] if they touch the edge 
            // of the domain in preparation for them being shrunk by n_error_buf[levf][idir] later.
            // We want the net effect to be that grids are NOT shrunk away
            // from the edges of the domain.
            //
            for (BoxList::iterator blt = bl_tagged.begin(), End = bl_tagged.end();
                 blt != End;
                 ++blt)
            {
                for (int idir = 0; idir < AMREX_SPACEDIM; idir++)
                {
                    if (blt->smallEnd(idir) == Geom(levf).Domain().smallEnd(idir))
                        blt->growLo(idir,n_error_buf[levf][idir]);
                    if (blt->bigEnd(idir) == Geom(levf).Domain().bigEnd(idir))
                        blt->growHi(idir,n_error_buf[levf][idir]);
                }
            }
            Box mboxF = amrex::grow(bl_tagged.minimalBox(),1);
            BoxList blFcomp;
            blFcomp.complementIn(mboxF,bl_tagged);
            blFcomp.simplify();
            bl_tagged.clear();

            const IntVect& iv = IntVect(AMREX_D_DECL(n_error_buf[levf][0]/ref_ratio[levf][0],
                                                     n_error_buf[levf][1]/ref_ratio[levf][1],
                                                     n_error_buf[levf][2]/ref_ratio[levf][2]));
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
            for (int idir = 0; idir < AMREX_SPACEDIM; idir++)
            {
                if (              n_error_buf[levf][idir] >  n_error_buf[levc][idir]*ref_ratio[levc][idir])
                    baF.grow(idir,n_error_buf[levf][idir]  - n_error_buf[levc][idir]*ref_ratio[levc][idir]);
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
        for (int n=0; n<AMREX_SPACEDIM; n++) {
            bl_max = std::max(bl_max,bf_lev[levc][n]);
        }
        if (bl_max >= 1) {
            tags.coarsen(bf_lev[levc]);
        } else {
            amrex::Abort("blocking factor is too small relative to ref_ratio");
        }
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
	Vector<IntVect> tagvec;
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
            if (use_new_chop)
            {
               clist.new_chop(grid_eff);
            } else {
               clist.chop(grid_eff);
            }
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
		    new_bx = amrex::intersect(new_bx,Geom(levc).Domain());
		}
	    }

            const IntVect& largest_grid_size = max_grid_size[levf] / ref_ratio[levc];
            //
            // Ensure new grid boxes are at most max_grid_size in index dirs.
            //
            new_bx.maxSize(largest_grid_size);

            //
            // Refine up to levf.
            //
            new_bx.refine(ref_ratio[levc]);
            BL_ASSERT(new_bx.isDisjoint());

	    if (new_bx.size()>0) {
		if ( !(Geom(levf).Domain().contains(BoxArray(new_bx).minimalBox())) ) {
		    new_bx = amrex::intersect(new_bx,Geom(levf).Domain());
		}
	    }

            if(levf > useFixedUpToLevel()) {
              new_grids[levf].define(new_bx);
	    }
        }
    }

    for (int lev = lbase+1; lev <= new_finest; ++lev) {
        if (new_grids[lev].empty())
        {
            if (!(useFixedCoarseGrids() && lev<useFixedUpToLevel()) ) {
                amrex::Abort("AmrMesh::MakeNewGrids: how did this happen?");
            }
        }
        else if (refine_grid_layout)
        {
            ChopGrids(lev,new_grids[lev],ParallelDescriptor::NProcs());
            if (new_grids[lev] == grids[lev]) {
                new_grids[lev] = grids[lev]; // to avoid dupliates
            }
        }
    }
}

void
AmrMesh::MakeNewGrids (Real time)
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
	Vector<BoxArray> new_grids(max_level+1);
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
        if (iterate_on_new_grids)
	{
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
}

void
AmrMesh::ProjPeriodic (BoxList& blout, const Geometry& geom)
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
    AMREX_D_TERM( nist , =njst , =nkst ) = -1;
    AMREX_D_TERM( niend , =njend , =nkend ) = +1;

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

void
AmrMesh::checkInput ()
{
    if (max_level < 0)
        amrex::Error("checkInput: max_level not set");

    //
    // Check level dependent values.
    //
    for (int i = 0; i < max_level; i++)
    {
        if (MaxRefRatio(i) < 2 || MaxRefRatio(i) > 12)
            amrex::Error("Amr::checkInput: bad ref_ratios");
    }

    const Box& domain = Geom(0).Domain();
    if (!domain.ok())
        amrex::Error("level 0 domain bad or not set");

    //
    // Check that domain size is a multiple of blocking_factor[0].
    //   (only check if blocking_factor <= max_grid_size)
    //
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
    {
        int len = domain.length(idim);
        if (blocking_factor[0][idim] <= max_grid_size[0][idim])
           if (len%blocking_factor[0][idim] != 0)
           {
              amrex::Print() << "domain size in direction " << idim << " is " << len << std::endl;
              amrex::Print() << "blocking_factor is " << blocking_factor[0][idim] << std::endl;
              amrex::Error("domain size not divisible by blocking_factor");
           }
    }

    //
    // Check that max_grid_size is a multiple of blocking_factor at every level.
    //   (only check if blocking_factor <= max_grid_size)
    //
    for (int i = 0; i < max_level; i++)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
           if (blocking_factor[i][idim] <= max_grid_size[i][idim])
              if (max_grid_size[i][idim]%blocking_factor[i][idim] != 0) {
              {
                 amrex::Print() << "max_grid_size in direction " << idim 
                                << " is " << max_grid_size[i][idim] << std::endl;
                 amrex::Print() << "blocking_factor is " << blocking_factor[i][idim] << std::endl;
                 amrex::Error("max_grid_size not divisible by blocking_factor");
              }
            }
        }
    }

    if( ! (Geometry::ProbDomain().volume() > 0.0) ) {
        amrex::Error("Amr::checkInput: bad physical problem size");
    }

    if(verbose > 0) {
        amrex::Print() << "Successfully read inputs file ... " << '\n';
    }
}

long
AmrMesh::CountCells (int lev) noexcept
{
    return grids[lev].numPts();
}


}
