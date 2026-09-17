// Gridding for odd refinement ratios and no_chop_dir. The legacy paths
// remain in AMReX_AmrMesh.cpp to preserve existing applications.
#include <AMReX.H>
#include <AMReX_AmrMesh.H>
#include <AMReX_Cluster.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#ifdef AMREX_USE_BITTREE
#include <AMReX_Bittree.H>
#endif
#include <optional>

namespace amrex {

IntVect
AmrMesh::effectiveMaxGridSize (int lev) const noexcept
{
    IntVect mgs = max_grid_size[lev];
    if (no_chop_dir >= 0) {
        AMREX_ASSERT(no_chop_dir < AMREX_SPACEDIM);
        // A grid can never be longer than the domain, so this makes the
        // max_grid_size constraint vacuous in that direction.
        mgs[no_chop_dir] = std::max(mgs[no_chop_dir],
                                    Geom(lev).Domain().length(no_chop_dir));
    }
    return mgs;
}

IntVect
AmrMesh::bfLev (int lev) const noexcept
{
    AMREX_ASSERT(lev >= 0 && lev < max_level);
    if (useLegacyGridding()) {
        return IntVect(1).max(blocking_factor[lev+1]/ref_ratio[lev]);
    }
    IntVect bf;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        int b = std::max(1, blocking_factor[lev+1][idim]/ref_ratio[lev][idim]);
        if (Geom(lev).isPeriodic(idim)) {
            // The periodic mapping of the coarsened tags requires the
            // domain to be divisible by the coarsening factor.
            int const len = Geom(lev).Domain().length(idim);
            if (len % b != 0) {
                int p = 1;
                while (2*p <= b && len % (2*p) == 0) { p *= 2; }
                b = p;
            }
        }
        bf[idim] = b;
    }
    return bf;
}

namespace {

// A "partial" box touches the upper domain boundary in idim, is longer
// than chunk there, and its length is not a multiple of unit.  This only
// happens when the domain is not divisible by the blocking factor.
bool isPartialBox (Box const& b, int idim, int chunk, int unit, Box const& domain)
{
    int const len = b.length(idim);
    return (len % unit != 0) && (len > chunk) && (b.bigEnd(idim) == domain.bigEnd(idim));
}

bool hasPartialBox (BoxArray const& ba, int idim, int chunk, int unit, Box const& domain)
{
    for (int i = 0; i < ba.size(); ++i) {
        if (isPartialBox(ba[i], idim, chunk, unit, domain)) { return true; }
    }
    return false;
}

// Chop the boxes in direction idim only, so that they are no longer than
// chunk there, keeping them coarsenable by rr.  A partial box (see above)
// is extended to a multiple of unit first, chopped, and clipped to the
// domain again, so that its pieces are multiples of unit with the
// remainder attached to the last one.  chunk must not chop the other
// directions.
void chopBoxes (BoxArray& ba, int idim, int rr, IntVect const& chunk, int unit,
                Box const& domain)
{
    IntVect crr(1);
    crr[idim] = rr;
    IntVect const cchunk = chunk / crr;
    BoxList bl(ba.ixType());
    for (int i = 0; i < ba.size(); ++i) {
        Box const& b = ba[i];
        int const len = b.length(idim);
        if (isPartialBox(b, idim, chunk[idim], unit, domain)) {
            int const rem = len % unit;
            Box be = b;
            be.growHi(idim, unit-rem);
            BoxList pieces(be);
            pieces.coarsen(crr);
            pieces.maxSize(cchunk);
            pieces.refine(crr);
            Vector<Box> v;
            for (auto const& p : pieces) {
                Box const q = p & domain;
                if (q.ok()) { v.push_back(q); }
            }
            if (v.size() >= 2 && v.back().length(idim) < unit) {
                v[v.size()-2].setBig(idim, v.back().bigEnd(idim));
                v.pop_back();
            }
            bl.join(v);
        } else {
            BoxList tmp(b);
            tmp.coarsen(crr);
            tmp.maxSize(cchunk);
            tmp.refine(crr);
            bl.join(tmp);
        }
    }
    ba.repartition(std::move(bl));
}

}

void
AmrMesh::ChopGridsExtended (int lev, BoxArray& ba, int target_size) const
{
    IntVect chop_dims = refine_grid_layout_dims;
    if (no_chop_dir >= 0) { chop_dims[no_chop_dir] = 0; }
    if (chop_dims == 0) { return; }

    // Use the effective max_grid_size here because ba is allowed to violate
    // max_grid_size in the no-chop direction, and the maxSize calls below
    // apply chunk in every direction, not just the one being refined.
    IntVect chunk = effectiveMaxGridSize(lev);
    chunk.min(Geom(lev).Domain().length());

    // Note that ba already satisfies the effective max_grid_size requirement
    // and it's coarsenable if it's a fine level BoxArray.

    while (ba.size() < target_size)
    {
        IntVect chunk_prev = chunk;

        std::array<std::pair<int,int>,AMREX_SPACEDIM>
            chunk_dir{AMREX_D_DECL(std::make_pair(chunk[0],int(0)),
                                   std::make_pair(chunk[1],int(1)),
                                   std::make_pair(chunk[2],int(2)))};
        std::ranges::sort(chunk_dir);

        for (int idx = AMREX_SPACEDIM-1; idx >= 0; idx--) {
            int idim = chunk_dir[idx].second;
            if (chop_dims[idim]) {
                int new_chunk_size = chunk[idim] / 2;
                int rr = (lev > 0) ? ref_ratio[lev-1][idim] : 1;
                if (rr > 1) {
                    new_chunk_size = (new_chunk_size/rr) * rr;
                }
                if (new_chunk_size != 0 &&
                    new_chunk_size%blocking_factor[lev][idim] == 0)
                {
                    chunk[idim] = new_chunk_size;
                    int const unit = (lev > 0) ? bfLev(lev-1)[idim]*rr : 1;
                    Box const& domain = Geom(lev).Domain();
                    // Chop in idim only.  A box at a truncated domain
                    // boundary may be longer than chunk in another
                    // direction, and must not be split there.
                    IntVect chunk1 = domain.length();
                    chunk1[idim] = chunk[idim];
                    bool partial_possible = false;
                    if (lev > 0 && domain.length(idim) % unit != 0) {
                        for (auto const& r : ref_ratio[lev-1]) {
                            if (r != 1 && (r%2 != 0)) { partial_possible = true; }
                        }
                    }
                    if (partial_possible && hasPartialBox(ba, idim, chunk[idim], unit, domain)) {
                        chopBoxes(ba, idim, rr, chunk1, unit, domain);
                    } else if (rr == 1) {
                        ba.maxSize(chunk1);
                    } else {
                        IntVect bf(1);
                        bf[idim] = rr;
                        ba.minmaxSize(bf, chunk1);
                    }
                    break;
                }
            }
        }

        if (chunk == chunk_prev) {
            break;
        }
    }
}

BoxArray
AmrMesh::MakeBaseGridsNoChop () const
{
    const Box& dom = geom[0].Domain();
    BoxArray ba;

    // With no_chop_dir and blocking factor 1 in the other directions, the
    // domain is decomposed into nearly equal pieces without any alignment
    // requirement, which works for any number of cells.
    bool use_decompose = (no_chop_dir >= 0);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (idim != no_chop_dir && blocking_factor[0][idim] != 1) {
            use_decompose = false;
        }
    }

    if (use_decompose)
    {
        Array<bool,AMREX_SPACEDIM> decomp;
        int nboxes = 1;
        bool chop_for_procs = false;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            decomp[idim] = (idim != no_chop_dir);
            if (decomp[idim]) {
                int const mgs = max_grid_size[0][idim];
                nboxes *= (dom.length(idim) + mgs - 1) / mgs;
                if (refine_grid_layout_dims[idim]) { chop_for_procs = true; }
            }
        }
        if (refine_grid_layout && chop_for_procs && ParallelDescriptor::NProcs() > nboxes) {
            nboxes = ParallelDescriptor::NProcs();
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                decomp[idim] = decomp[idim] && refine_grid_layout_dims[idim];
            }
        }
        ba = amrex::decompose(dom, nboxes, decomp);
        // decompose only honors the total count; enforce max_grid_size
        // direction by direction.
        ba.maxSize(effectiveMaxGridSize(0));
    }
    else
    {
        IntVect fac(2);
        const Box dom2 = amrex::refine(amrex::coarsen(dom,2),2);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (dom.length(idim) != dom2.length(idim)) {
                fac[idim] = 1;
            }
        }
        ba = BoxArray(amrex::coarsen(dom,fac));
        ba.maxSize(effectiveMaxGridSize(0)/fac);
        ba.refine(fac);
        // Boxes in ba have even number of cells in each direction
        // unless the domain has odd number of cells in that direction.
        if (refine_grid_layout) {
            ChopGrids(0, ba, ParallelDescriptor::NProcs());
        }
    }

    if (ba == grids[0]) {
        ba = grids[0];  // to avoid duplicates
    }
    PostProcessBaseGrids(ba);
    return ba;
}

void
AmrMesh::MakeNewGridsExtended (int lbase, Real time, int& new_finest, Vector<BoxArray>& new_grids)
{
    BL_PROFILE("AmrMesh::MakeNewGrids()");

    BL_ASSERT(lbase < max_level);

    // Add at most one new level
    int max_crse = std::min(finest_level, max_level-1);

    if (new_grids.size() < max_crse+2) { new_grids.resize(max_crse+2); }

#ifdef AMREX_USE_BITTREE
    if(!use_bittree) {
#endif

    //
    // Construct problem domain at each level.
    //
    Vector<IntVect> bf_lev(max_level); // Blocking factor at each level.
    Vector<IntVect> rr_lev(max_level);
    Vector<Box>     pc_domain(max_level);  // Coarsened problem domain.

    for (int i = 0; i <= max_crse; i++)
    {
        bf_lev[i] = bfLev(i);
    }
    for (int i = lbase; i < max_crse; i++)
    {
        for (int n=0; n<AMREX_SPACEDIM; n++) {
            // For odd ratios checkInput makes sure this is an integer; even
            // ratios keep the old integer division.
            rr_lev[i][n] = (ref_ratio[i][n]*bf_lev[i][n])/bf_lev[i+1][n];
        }
    }
    for (int i = lbase; i <= max_crse; i++) {
        pc_domain[i] = amrex::coarsen(Geom(i).Domain(),bf_lev[i]);
    }
    //
    // Construct proper nesting domains.
    //
    Vector<BoxArray> p_n_ba(max_level); // Proper nesting domain.
    Vector<BoxArray> p_n_comp_ba(max_level); // Complement proper nesting domain.
    BoxList p_n, p_n_comp;

    BoxList bl_base = grids[lbase].simplified_list();
    bl_base.coarsen(bf_lev[lbase]);
    p_n_comp.parallelComplementIn(pc_domain[lbase],bl_base);
    bl_base.clear();
    p_n_comp.simplify();
    p_n_comp.accrete(n_proper);
    if (geom[lbase].isAnyPeriodic()) {
        ProjPeriodic(p_n_comp, pc_domain[lbase], geom[lbase].isPeriodic());
    }

    p_n_comp_ba[lbase].define(std::move(p_n_comp));
    p_n_comp = BoxList();

    p_n.parallelComplementIn(pc_domain[lbase],p_n_comp_ba[lbase]);
    p_n.simplify();

    p_n_ba[lbase].define(std::move(p_n));
    p_n = BoxList();

    for (int i = lbase+1; i <= max_crse; i++)
    {
        p_n_comp = p_n_comp_ba[i-1].boxList();

        // Need to simplify p_n_comp or the number of grids can too large for many levels.
        p_n_comp.simplify();

        p_n_comp.refine(rr_lev[i-1]);
        p_n_comp.accrete(n_proper);

        if (geom[i].isAnyPeriodic()) {
            ProjPeriodic(p_n_comp, pc_domain[i], geom[i].isPeriodic());
        }

        p_n_comp_ba[i].define(std::move(p_n_comp));
        p_n_comp = BoxList();

        p_n.parallelComplementIn(pc_domain[i],p_n_comp_ba[i]);
        p_n.simplify();

        p_n_ba[i].define(std::move(p_n));
        p_n = BoxList();
    }

    //
    // Now generate grids from finest level down.
    //
    new_finest = lbase;

    for (int levc = max_crse; levc >= lbase; levc--)
    {
        int levf = levc+1;

        TagBoxArray tags(grids[levc],dmap[levc],n_error_buf[levc]);

        //
        // Only use error estimation to tag cells for the creation of new grids
        //      if the grids at that level aren't already fixed.
        //

        if ( ! (useFixedCoarseGrids() && levc < useFixedUpToLevel()) ) {
            ErrorEst(levc, tags, time, 0);
        }

        //
        // Buffer error cells.
        //
        tags.buffer(n_error_buf[levc]);

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
            // Fixed or caller-supplied grids may have arbitrary cuts.
            bool const may_overlap = !grids[levc].coarsenable(bf_lev[levc]);
            tags.coarsenMayOverlap(bf_lev[levc], may_overlap);
        } else {
            amrex::Abort("blocking factor is too small relative to ref_ratio");
        }
        //
        // Remove or add tagged points which violate/satisfy additional
        // user-specified criteria.
        //
        ManualTagsPlacement(levc, tags, bf_lev);
        //
        // If new grids have been constructed above this level, project
        // those grids down and tag cells on intersections to ensure proper
        // nesting. Note that the projected BoxArray may contain cells
        // outside the TagBoxArray's domain. For those, we collect them in
        // Vector tag_proj.
        //
        Vector<IntVect> tags_proj;
        if (levf < new_finest) {
            BoxArray ba_tags = tags.boxArray();
            BoxList bl_proj = new_grids[levf+1].simplified_list();
            auto& bxs = bl_proj.data();
            Long nbxs = bl_proj.size();
            Vector<Vector<IntVect>> tags_proj_priv(OpenMP::get_max_threads());
            Box domain = Geom(levc).Domain();
            domain.coarsen(bf_lev[levc]);
            auto const& domain_lo = domain.smallEnd();
            auto const& domain_hi = domain.bigEnd();
            auto const& domain_len = domain.length();
            auto const& is_periodic = Geom(levc).isPeriodicArray();
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
            {
                auto& tv = tags_proj_priv[OpenMP::get_thread_num()];
                std::vector<std::pair<int,Box>> isects;
                BoxList bl1(ba_tags.ixType());
                BoxList bl2(ba_tags.ixType());
                BoxList bltmp;
#ifdef AMREX_USE_OMP
#pragma omp for
#endif
                for (Long ibx = 0; ibx < nbxs; ++ibx) {
                    Box& b = bxs[ibx];
                    b.coarsen(ref_ratio[levf]);
                    b.grow(bf_lev[levf]*n_proper);
                    b.coarsen(ref_ratio[levc]);
                    b.coarsen(bf_lev[levc]);

                    ba_tags.intersections(b, isects, false, tags.nGrowVect());
                    bl1.clear();
                    bl1.push_back(b);
                    for (auto const& kv : isects) {
                        bl2.clear();
                        for (auto & btmp : bl1) {
                            amrex::boxDiff(bltmp, btmp, kv.second);
                            bl2.join(bltmp);
                        }
                        std::swap(bl1,bl2);
                    }
                    for (auto const& bleft : bl1) {
                        amrex::LoopOnCpu(bleft, [&] (IntVect iv)
                        {
                            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                                if (iv[idim] < domain_lo[idim]) {
                                    if (is_periodic[idim]) {
                                        iv[idim] += domain_len[idim];
                                    } else {
                                        return;
                                    }
                                } else if (iv[idim] > domain_hi[idim]) {
                                    if (is_periodic[idim]) {
                                        iv[idim] -= domain_len[idim];
                                    } else {
                                        return;
                                    }
                                }
                            }
                            tv.push_back(iv);
                        });
                    }
                }
            }

            BoxArray ba_proj(std::move(bl_proj));
            tags.setVal(ba_proj,TagBox::SET);

            Long nextra = 0;
            for (auto const& tv : tags_proj_priv) {
                nextra += tv.size();
            }
            if (nextra > 0) {
                tags_proj.reserve(nextra);
                for (auto const& tv : tags_proj_priv) {
                    tags_proj.insert(std::end(tags_proj), std::begin(tv), std::end(tv));
                }
                amrex::RemoveDuplicates(tags_proj);
            }
        }
        //
        // Map tagged points through periodic boundaries, if any.
        //
        tags.mapPeriodicRemoveDuplicates(Geometry(pc_domain[levc],
                                                  Geom(levc).ProbDomain(),
                                                  Geom(levc).CoordInt(),
                                                  Geom(levc).isPeriodic()));
        //
        // Remove cells outside proper nesting domain for this level.
        //
        tags.setVal(p_n_comp_ba[levc],TagBox::CLEAR);
        p_n_comp_ba[levc].clear();
        //
        // Create initial cluster containing all tagged points.
        //
        Gpu::PinnedVector<IntVect> tagvec;
        tags.collate(tagvec);
        tags.clear();

        if (!tags_proj.empty()) {
            tagvec.insert(tagvec.end(), tags_proj.data(), tags_proj.data()+tags_proj.size());
            tags_proj.clear();
        }

        if (!tagvec.empty())
        {
            //
            // Created new level, now generate efficient grids.
            //
            if ( !(useFixedCoarseGrids() && levc<useFixedUpToLevel()) ) {
                new_finest = std::max(new_finest,levf);
            }

            if (levf > useFixedUpToLevel()) {
                bool odd_ref_ratio = false;
                for (auto const& rr : ref_ratio[levc]) {
                    if (rr != 1 && (rr%2 != 0)) {
                        odd_ref_ratio = true;
                    }
                }

                BoxList new_bx;
                if (ParallelDescriptor::IOProcessor()) {
                    BL_PROFILE("AmrMesh-cluster");
                    //
                    // Construct initial cluster.
                    //
                    ClusterList clist(tagvec.data(), static_cast<Long>(tagvec.size()),
                                      pc_domain[levc], refine_whole_domain_dir);
                    if (use_new_chop) {
                        clist.new_chop(grid_eff);
                    } else {
                        clist.chop(grid_eff);
                    }
                    // Only odd ratios need a copy for the boundary extension.
                    std::optional<BoxArray> p_n_copy;
                    if (odd_ref_ratio) { p_n_copy = p_n_ba[levc]; }
                    clist.intersect(p_n_ba[levc]);
                    //
                    // Efficient properly nested Clusters have been constructed
                    // now generate list of grids at level levf.
                    //
                    clist.boxList(new_bx);
                    new_bx.simplify();

                    // Pair a partial last cell with its neighbor before removing
                    // overlaps, so the extension cannot be trimmed back to one cell.
                    if (odd_ref_ratio)
                    {
                        Box const& pcd = pc_domain[levc];
                        IntVect paired(0);
                        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                            if (Geom(levc).Domain().length(idim) % bf_lev[levc][idim] != 0
                                && pcd.length(idim) > 1)
                            {
                                paired[idim] = 1;
                            }
                        }
                        bool thin = false;
                        for (auto const& b : new_bx) {
                            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                                if (paired[idim] && b.bigEnd(idim) == pcd.bigEnd(idim)
                                    && b.length(idim) == 1) { thin = true; }
                            }
                        }
                        if (thin) {
                            // Temporarily identify the last two cells. Disjoint boxes
                            // in this space stay disjoint when the pair is expanded.
                            for (auto& b : new_bx) {
                                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                                    if (paired[idim]) {
                                        int const hi = pcd.bigEnd(idim)-1;
                                        b.setSmall(idim, std::min(b.smallEnd(idim), hi));
                                        b.setBig(idim, std::min(b.bigEnd(idim), hi));
                                    }
                                }
                            }
                            BoxArray ba(std::move(new_bx));
                            ba.removeOverlap(false);
                            new_bx = ba.boxList();
                            BoxList nested(new_bx.ixType());
                            for (auto b : new_bx) {
                                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                                    if (paired[idim] && b.bigEnd(idim) == pcd.bigEnd(idim)-1) {
                                        b.setBig(idim, pcd.bigEnd(idim));
                                    }
                                }
                                if (p_n_copy->contains(b)) {
                                    nested.push_back(b);
                                } else {
                                    nested.join(amrex::intersect(*p_n_copy, b).boxList());
                                }
                            }
                            new_bx = std::move(nested);
                            new_bx.simplify();
                        }
                    }

                    if (no_chop_dir >= 0) {
                        // No two grids may share a face normal to no_chop_dir.
                        // Nothing after this point chops in that direction.
                        new_bx.mergeAlongDir(no_chop_dir);
                    }
                }
                new_bx.Bcast();  // Broadcast the new BoxList to other processes

                // The boxes are in the index space of level levc coarsened
                // by bf_lev[levc].

                if (odd_ref_ratio)
                {
                    // This approach imposes max_grid_size (suitably scaled) before
                    //     refining so as to ensure fine grids align with coarse grids

                    // In each direction where max_grid_size is a multiple
                    // of the grid unit bf_lev*ref_ratio, chop before
                    // refining by bf_lev so that the grids are multiples of
                    // the unit even when the domain is not divisible by it.
                    // Other directions are chopped after refining, as we
                    // have always done.  The boxes are inside pc_domain, so
                    // using its length as chunk means no chop.
                    IntVect const emgs = effectiveMaxGridSize(levf);
                    IntVect chunk = pc_domain[levc].length();
                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                        int bf = bf_lev[levc][idim];
                        int unit = bf * ref_ratio[levc][idim];
                        if (idim != no_chop_dir && ((bf & (bf-1)) == 0)
                            && (emgs[idim]%unit == 0))
                        {
                            chunk[idim] = emgs[idim] / unit;
                        }
                    }
                    // Chop as before, except that a box at the upper
                    // boundary of a direction with a partial last cell is
                    // split with the larger pieces last if the usual split
                    // would leave that cell alone as a thin grid.
                    {
                        Box const& pcd = pc_domain[levc];
                        BoxList chopped(new_bx.ixType());
                        Vector<Box> pieces, next;
                        for (auto const& b : new_bx) {
                            BoxList plain(b);
                            plain.maxSize(chunk);
                            IntVect fix(0);
                            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                                if (Geom(levc).Domain().length(idim) % bf_lev[levc][idim] != 0
                                    && b.bigEnd(idim) == pcd.bigEnd(idim)
                                    && b.length(idim) > chunk[idim])
                                {
                                    for (auto const& q : plain) {
                                        if (q.bigEnd(idim) == b.bigEnd(idim) && q.length(idim) == 1) {
                                            fix[idim] = 1;
                                        }
                                    }
                                }
                            }
                            if (fix == 0) {
                                chopped.join(plain);
                                continue;
                            }
                            IntVect c = chunk;
                            pieces.assign(1, b);
                            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                                if (!fix[idim]) { continue; }
                                next.clear();
                                for (auto const& p : pieces) {
                                    int const len = p.length(idim);
                                    int const n = (len + chunk[idim] - 1) / chunk[idim];
                                    int const base = len / n;
                                    int const extra = len - base*n;
                                    int lo = p.smallEnd(idim);
                                    for (int k = 0; k < n; ++k) {
                                        int const l = base + ((k >= n-extra) ? 1 : 0);
                                        Box q = p;
                                        q.setSmall(idim, lo);
                                        q.setBig(idim, lo+l-1);
                                        lo += l;
                                        next.push_back(q);
                                    }
                                }
                                pieces.swap(next);
                                c[idim] = b.length(idim); // done in this direction
                            }
                            BoxList bl(new_bx.ixType());
                            bl.join(pieces);
                            bl.maxSize(c);
                            chopped.join(bl);
                        }
                        new_bx = std::move(chopped);
                    }

                    new_bx.refine(bf_lev[levc]);
                    if (new_bx.size()>0) {
                        // Chop new grids outside domain
                        new_bx.intersect(Geom(levc).Domain());
                    }

                    //
                    // Impose max_grid_size (suitably coarsened)
                    //
                    AMREX_ASSERT(emgs.allGE(ref_ratio[levc]));
                    new_grids[levf] = BoxArray(std::move(new_bx),
                                               effectiveMaxGridSize(levf)/ref_ratio[levc]);

                    //
                    // Refine up to levf.
                    //
                    new_grids[levf].refine(ref_ratio[levc]);
                }
                else
                {
                    // This approach imposes max_grid_size after refining.
                    // For ref_ratio = 3 this can create fine grids that do not correctly divide by 3,
                    //     but we leave it here so as not to change the gridding in
                    //     existing ref_ratio = 2 or 4 applications

                    new_bx.refine(bf_lev[levc]);
                    if (new_bx.size()>0) {
                        // Chop new grids outside domain
                        new_bx.intersect(Geom(levc).Domain());
                    }

                    //
                    // Refine up to levf.
                    //
                    new_bx.refine(ref_ratio[levc]);

                    //
                    // Impose max_grid_size
                    //
                    new_grids[levf] = BoxArray(std::move(new_bx), effectiveMaxGridSize(levf));
                }
                BL_ASSERT(new_grids[levf].isDisjoint());
            }
        }
    }

#if 0
    if (!useFixedCoarseGrids()) {
        // check proper nesting
        // This check does not consider periodic boundary and could fail if
        // the blocking factor is not the same on all levels.
        for (int lev = lbase+1; lev <= new_finest; ++lev) {
            BoxArray const& cba = (lev == lbase+1) ? grids[lev-1] : new_grids[lev-1];
            BoxArray const& fba = amrex::coarsen(new_grids[lev],ref_ratio[lev-1]);
            IntVect np = bf_lev[lev-1] * n_proper;
            Box const& cdomain = Geom(lev-1).Domain();
            for (int i = 0, N = fba.size(); i < N; ++i) {
                Box const& fb = amrex::grow(fba[i],np) & cdomain;
                if (!cba.contains(fb,true)) {
                    amrex::Abort("AmrMesh::MakeNewGrids: new grids not properly nested");
                }
            }
        }
    }
#endif

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
                new_grids[lev] = grids[lev]; // to avoid duplicates
            }
        }
    }

#ifdef AMREX_USE_BITTREE
    }
#endif

#ifdef AMREX_USE_BITTREE
    // Bittree version
    if(use_bittree) {
        // Initialize BT refinement
        btmesh->refine_init();

        // -------------------------------------------------------------------
        // Use tagging data to mark BT for refinement, then use the new bitmap
        // to calculate the new grids.
        auto tree0 = btmesh->getTree();

        // [1] Error Estimation and tagging
        // btTags is indexed by bitid, Bittree's internal indexing scheme.
        // For any id, btTags = 1 if should be parent, -1 if should not be parent (or not exist).
        std::vector<int> btTags(tree0->id_upper_bound(),0);

        for (int lev=max_crse; lev>=lbase; --lev) {

            TagBoxArray tags(grids[lev],dmap[lev], n_error_buf[lev]);
            ErrorEst(lev, tags, time, 0);
            tags.buffer(n_error_buf[lev]);

            for (MFIter mfi(tags); mfi.isValid(); ++mfi) {
                auto const& tagbox = tags.const_array(mfi);
                bool has_set_tags = amrex::Reduce::AnyOf(mfi.validbox(),
                                                         [=] AMREX_GPU_DEVICE (int i, int j, int k)
                                                         {
                                                              return tagbox(i,j,k)!=TagBox::CLEAR;
                                                         });

                // Set the values of btTags.
                int bitid = btUnit::getBitid(btmesh.get(),false,lev,mfi.index());
                // TODO Check lev == tree0->block_level(bitid)
                if(has_set_tags) {
                    btTags[bitid] = 1;
                }
                else {
                    btTags[bitid] = -1;
                }
            }
        }

        // [2] btRefine - check for proper octree nesting and update bitmap
        MPI_Comm comm = ParallelContext::CommunicatorSub();
        int changed = btUnit::btRefine(btmesh.get(), btTags, max_crse, lbase, grids, dmap, comm);

        // [3] btCalculateGrids - use new bitmap to generate new grids
        if (changed>0) {
            btUnit::btCalculateGrids(btmesh.get(),lbase,new_finest,new_grids,max_grid_size);
        } else {
            new_finest = finest_level;
            for(int i=0; i<=finest_level; ++i) {
                new_grids[i] = grids[i];
            }
        }

        // Finalize BT refinement
        btmesh->refine_apply();
    }
#endif

}

void
AmrMesh::checkInputExtended ()
{
    if (max_level < 0) {
        amrex::Error("checkInput: max_level not set");
    }

    //
    // Check level dependent values.
    //
    for (int i = 0; i < max_level; i++)
    {
        if (MaxRefRatio(i) < 2) {
            amrex::Warning("Amr::checkInput: ref_ratios all equal to one!");
        }
    }

    const Box& domain = Geom(0).Domain();
    if (!domain.ok()) {
        amrex::Error("level 0 domain bad or not set");
    }

    //
    // Check the direction in which the grids are never chopped.
    //
    if (no_chop_dir >= 0)
    {
        if (no_chop_dir >= AMREX_SPACEDIM) {
            amrex::Error("Amr::checkInput: no_chop_dir is out of range");
        }
#ifdef AMREX_USE_BITTREE
        if (use_bittree) {
            amrex::Error("Amr::checkInput: no_chop_dir does not work with bittree");
        }
#endif
    }

    //
    // Check that domain size is a multiple of blocking_factor[0].
    //   (only check if blocking_factor <= max_grid_size, and not in
    //   no_chop_dir where blocking_factor and max_grid_size are ignored)
    //
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
    {
        int len = domain.length(idim);
        if (idim != no_chop_dir && blocking_factor[0][idim] <= max_grid_size[0][idim]) {
            if (len%blocking_factor[0][idim] != 0)
            {
                amrex::Print() << "domain size in direction " << idim << " is " << len << '\n'
                               << "blocking_factor is " << blocking_factor[0][idim] << '\n';
                amrex::Error("domain size not divisible by blocking_factor");
            }
        }
    }

    auto is_pow2 = [] (int k) { return k > 0 && (k & (k-1)) == 0; };

    // Grids on level i > 0 are built from tags on level i-1 coarsened by
    // bfLev(i-1) and refined back, so they are multiples of the "grid unit"
    // bfLev(i-1)*ref_ratio[i-1].
    auto grid_unit = [&] (int i, int idim) {
        return bfLev(i-1)[idim] * ref_ratio[i-1][idim];
    };

    //
    // Check that blocking_factor is a power of 2.  With an odd refinement
    // ratio, the blocking factor on a fine level may also be the ref ratio
    // times a power of 2 (e.g., 24 for ref_ratio 3).
    //
    for (int i = 0; i <= max_level; i++)
    {
        bool odd_rr = false;
        if (i > 0) {
            for (auto const& rr : ref_ratio[i-1]) {
                if (rr != 1 && (rr%2 != 0)) { odd_rr = true; }
            }
        }
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            int bf = blocking_factor[i][idim];
            bool ok = is_pow2(bf);
            if (!ok && odd_rr) {
                int rr = ref_ratio[i-1][idim];
                ok = (bf%rr == 0) && is_pow2(bf/rr);
            }
            if (!ok) {
                amrex::Print() << "blocking_factor on level " << i << " in direction "
                               << idim << " is " << bf << '\n';
                amrex::Error("Amr::checkInput: blocking_factor not power of 2 (or, for odd ref_ratio, ref_ratio times power of 2). You can bypass this by setting ParmParse runtime parameter amr.check_input=0, although we do not recommend it.");
            }
        }
    }

    //
    // Warn if the blocking factor cannot be satisfied with the ref ratio.
    //
    for (int i = 1; i <= max_level; i++) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            int bf = blocking_factor[i][idim];
            int rr = ref_ratio[i-1][idim];
            int unit = std::max(1, bf/rr) * rr; // before any periodic reduction
            if (unit%bf != 0) {
                amrex::Print() << "WARNING: blocking_factor " << bf << " on level " << i
                               << " in direction " << idim << " cannot be satisfied with ref_ratio "
                               << ref_ratio[i-1][idim] << ". Grids will be multiples of "
                               << unit << " instead.\n";
            }
        }
    }

    //
    // Check that max_grid_size is a multiple of blocking_factor at every level.
    //   (only check if max_level > 0 && blocking_factor <= max_grid_size)
    //
    if (max_level > 0) {
        for (int i = 0; i <= max_level; i++)
        {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if (idim != no_chop_dir && blocking_factor[i][idim] <= max_grid_size[i][idim]) {
                    if (max_grid_size[i][idim]%blocking_factor[i][idim] != 0) {
                        amrex::Print() << "max_grid_size in direction " << idim
                                       << " is " << max_grid_size[i][idim] << '\n'
                                       << "blocking_factor is " << blocking_factor[i][idim] << '\n';
                        amrex::Error("max_grid_size not divisible by blocking_factor");
                    }
                }
            }
        }
    }

    //
    // With a domain that is not divisible by the blocking factor, the grid
    // at the upper boundary can only be kept at least as thick as the
    // blocking factor if max_grid_size allows grids of two blocking factors.
    // This is only done for odd refinement ratios; even ratios keep the
    // old behavior.
    //
    for (int i = 1; i <= max_level; i++) {
        bool odd_rr = false;
        for (auto const& rr : ref_ratio[i-1]) {
            if (rr != 1 && (rr%2 != 0)) { odd_rr = true; }
        }
        if (!odd_rr) { continue; }
        IntVect const emgs = effectiveMaxGridSize(i);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            int const bf_lev = bfLev(i-1)[idim];
            int const unit = grid_unit(i,idim);
            if (idim != no_chop_dir && Geom(i-1).Domain().length(idim) % bf_lev != 0
                && emgs[idim] < 2*unit)
            {
                amrex::Print() << "On level " << i << " in direction " << idim
                               << " max_grid_size is " << emgs[idim] << " and the grids are multiples of "
                               << unit << ", but the level " << i-1 << " domain size "
                               << Geom(i-1).Domain().length(idim) << " is not divisible by "
                               << bf_lev << ".\n";
                amrex::Error("max_grid_size must be at least twice the blocking factor when the domain is not divisible by it");
            }
        }
    }

    //
    // Check that blocking_factor does not vary too much between levels.
    // Grids on level i (1 <= i < max_level) must be coarsenable by
    // bf_lev[i] = max(1, blocking_factor[i+1]/ref_ratio[i]).  Level 0 is
    // exempt because its grids cover the whole domain.
    //
    for (int i = 1; i < max_level; i++) {
        IntVect const bf_lev = bfLev(i);
        IntVect const emgs = effectiveMaxGridSize(i);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            // For even ratios this is the old check.
            bool odd_prev = false;
            for (auto const& rr : ref_ratio[i-1]) {
                if (rr != 1 && (rr%2 != 0)) { odd_prev = true; }
            }
            int const gu = odd_prev ? grid_unit(i,idim) : blocking_factor[i][idim];
            int unit = std::min(gu, emgs[idim]);
            if (unit % bf_lev[idim] != 0) {
                amrex::Print() << "Blocking factors on levels " << i << " and " << i+1
                               << " are " << blocking_factor[i] << " " << blocking_factor[i+1]
                               << ". Ref ratio is " << ref_ratio[i]
                               << ".  They vary too much between levels." << '\n';
                amrex::Error("Blocking factors vary too much between levels");
            }
        }
    }

    //
    // In periodic directions, the tag coarsening factor is reduced if the
    // domain is not divisible by it.
    //
    for (int i = 0; i < max_level; ++i) {
        IntVect const bf_lev = bfLev(i);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            int nominal = std::max(1, blocking_factor[i+1][idim]/ref_ratio[i][idim]);
            if (bf_lev[idim] != nominal) {
                amrex::Print() << "WARNING: Level " << i << " domain size in periodic direction "
                               << idim << " is " << Geom(i).Domain().length(idim)
                               << ", which is not divisible by blocking_factor/ref_ratio = "
                               << nominal << ". Tags will be coarsened by " << bf_lev[idim]
                               << " instead, and level " << i+1 << " grids will be multiples of "
                               << bf_lev[idim]*ref_ratio[i][idim] << " in that direction.\n";
            }
        }
    }

    //
    // Check the direction in which the fine levels cover the entire domain.
    //
    if (refine_whole_domain_dir >= 0)
    {
        const int idim = refine_whole_domain_dir;
        if (idim >= AMREX_SPACEDIM) {
            amrex::Error("Amr::checkInput: refine_whole_domain_dir is out of range");
        }
#ifdef AMREX_USE_BITTREE
        if (use_bittree) {
            amrex::Error("Amr::checkInput: refine_whole_domain_dir does not work with bittree");
        }
#endif
    }

    if( ! (Geom(0).ProbDomain().volume() > 0.0) ) {
        amrex::Error("Amr::checkInput: bad physical problem size");
    }

    if(verbose > 0) {
        amrex::Print() << "Successfully read inputs file ... " << '\n';
    }
}

}
