// Exercise AmrMesh grid creation and check the properties of the resulting
// grids: max_grid_size, blocking factor alignment, coarsenability by
// ref_ratio, coverage of tagged cells, proper nesting, the no_chop_dir
// guarantees, and FillPatchTwoLevels on the resulting grids.

#include <AMReX.H>
#include <AMReX_AmrMesh.H>
#include <AMReX_BCRec.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_TagBox.H>
#include <AMReX_iMultiFab.H>

#include <algorithm>
#include <cmath>
#include <fstream>

using namespace amrex;

namespace {

class TestMesh : public AmrMesh
{
public:
    TestMesh ()
    {
        ParmParse pp("test");
        Vector<Real> c;
        pp.queryarr("centers", c);
        AMREX_ALWAYS_ASSERT(c.size() % AMREX_SPACEDIM == 0);
        for (int i = 0; i < c.size(); i += AMREX_SPACEDIM) {
            RealVect rv;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { rv[d] = c[i+d]; }
            m_centers.push_back(rv);
        }
        pp.queryarr("radius", m_radius);
        AMREX_ALWAYS_ASSERT(m_radius.size() >= 1);
    }

    // A cell is tagged if its center is within the radius of one of the
    // centers.  The radius shrinks on finer levels.
    [[nodiscard]] bool tagged (int lev, IntVect const& iv) const
    {
        Real const rad = m_radius[std::min(lev, int(m_radius.size())-1)];
        auto const& lgeom = Geom(lev);
        auto const dx = lgeom.CellSizeArray();
        auto const plo = lgeom.ProbLoArray();
        for (auto const& c : m_centers) {
            Real r2 = 0;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                Real x = plo[d] + (iv[d]+Real(0.5))*dx[d] - c[d];
                r2 += x*x;
            }
            if (r2 <= rad*rad) { return true; }
        }
        return false;
    }

    void ErrorEst (int lev, TagBoxArray& tags, Real /*time*/, int /*ngrow*/) override
    {
        for (MFIter mfi(tags); mfi.isValid(); ++mfi) {
            auto const& arr = tags.array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k)
            {
                if (tagged(lev, IntVect(AMREX_D_DECL(i,j,k)))) {
                    arr(i,j,k) = TagBox::SET;
                }
            });
        }
    }

    // Grids on level lev are multiples of this (except at the domain boundary).
    [[nodiscard]] IntVect gridUnit (int lev) const
    {
        IntVect g(1);
        if (lev > 0) {
            g = bfLev(lev-1) * refRatio(lev-1);
        }
        return g;
    }

    [[nodiscard]] int noChopDir () const { return no_chop_dir; }
    [[nodiscard]] bool refineGridLayout () const { return refine_grid_layout; }

private:
    Vector<RealVect> m_centers;
    Vector<Real> m_radius;
};

int nfail = 0;

void fail (std::string const& msg)
{
    amrex::Print() << "FAIL: " << msg << '\n';
    ++nfail;
}

// No two boxes share a face normal to dir.
void check_no_faces (BoxArray const& ba, int dir, int lev)
{
    std::vector<std::pair<int,Box>> isects;
    for (int i = 0; i < ba.size(); ++i) {
        Box s = ba[i];
        s.shift(dir, 1);
        ba.intersections(s, isects);
        for (auto const& is : isects) {
            if (is.first != i) {
                fail("level " + std::to_string(lev) + " boxes " + std::to_string(i)
                     + " and " + std::to_string(is.first) + " share a face normal to direction "
                     + std::to_string(dir));
            }
        }
    }
}

// The projections of the boxes along dir do not overlap.  This is what ERF
// requires for its vertical solves.
void check_not_split (BoxArray const& ba, int dir, int lev)
{
    BoxList bl(ba.ixType());
    for (int i = 0; i < ba.size(); ++i) {
        Box b = ba[i];
        b.setRange(dir, 0);
        bl.push_back(b);
    }
    if (!BoxArray(std::move(bl)).isDisjoint()) {
        fail("level " + std::to_string(lev) + " grids are split in direction "
             + std::to_string(dir));
    }
}

void check_level0 (TestMesh const& mesh)
{
    BoxArray const& ba = mesh.boxArray(0);
    Box const& domain = mesh.Geom(0).Domain();
    int const dir = mesh.noChopDir();

    if (!ba.isDisjoint()) { fail("level 0 grids overlap"); }
    if (ba.numPts() != domain.numPts() || !domain.contains(ba.minimalBox())) {
        fail("level 0 grids do not cover the domain");
    }
    IntVect const mgs = mesh.effectiveMaxGridSize(0);
    for (int i = 0; i < ba.size(); ++i) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (ba[i].length(d) > mgs[d]) {
                fail("level 0 box " + std::to_string(i)
                     + " longer than max_grid_size in direction " + std::to_string(d));
            }
        }
    }

    if (dir < 0) { return; }

    for (int i = 0; i < ba.size(); ++i) {
        if (ba[i].length(dir) != domain.length(dir)) {
            fail("level 0 box " + std::to_string(i) + " does not span the domain in direction "
                 + std::to_string(dir));
        }
    }
    check_not_split(ba, dir, 0);

    bool bf1 = true;
    int nboxes = 1;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (d != dir) {
            if (mesh.blockingFactor(0)[d] != 1) { bf1 = false; }
            int const mgs = mesh.maxGridSize(0)[d];
            nboxes *= (domain.length(d) + mgs - 1) / mgs;
        }
    }
    if (bf1 && mesh.refineGridLayout()) {
        nboxes = std::max(nboxes, ParallelDescriptor::NProcs());
    }
    if (bf1 && ba.size() < nboxes) {
        fail("level 0 has " + std::to_string(ba.size()) + " grids, expected at least "
             + std::to_string(nboxes));
    }
}

void check_level (TestMesh const& mesh, int lev)
{
    BoxArray const& ba = mesh.boxArray(lev);
    Box const& domain = mesh.Geom(lev).Domain();
    IntVect const rr = mesh.refRatio(lev-1);
    IntVect const unit = mesh.gridUnit(lev);
    IntVect const mgs = mesh.effectiveMaxGridSize(lev);

    amrex::Print() << "Level " << lev << ": " << ba.size() << " grids, "
                   << ba.numPts() << " cells, unit " << unit << '\n';

    if (!ba.isDisjoint()) { fail("level " + std::to_string(lev) + " grids overlap"); }
    if (!ba.coarsenable(rr)) {
        fail("level " + std::to_string(lev) + " grids not coarsenable by ref_ratio");
    }
    if (!domain.contains(ba.minimalBox())) {
        fail("level " + std::to_string(lev) + " grids outside domain");
    }

    for (int i = 0; i < ba.size(); ++i) {
        Box const& b = ba[i];
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (b.length(d) > mgs[d]) {
                fail("level " + std::to_string(lev) + " box " + std::to_string(i)
                     + " longer than max_grid_size in direction " + std::to_string(d));
            }
            // Boxes at a truncated upper boundary are extended inward, so
            // no box is thinner than the grid unit.
            if (b.length(d) < unit[d]) {
                fail("level " + std::to_string(lev) + " box " + std::to_string(i)
                     + " thinner than the blocking factor in direction " + std::to_string(d));
            }
            if ((b.smallEnd(d) - domain.smallEnd(d)) % unit[d] != 0) {
                fail("level " + std::to_string(lev) + " box " + std::to_string(i)
                     + " lower end not aligned to blocking factor in direction "
                     + std::to_string(d));
            }
            if ((b.bigEnd(d) + 1 - domain.smallEnd(d)) % unit[d] != 0
                && b.bigEnd(d) != domain.bigEnd(d))
            {
                fail("level " + std::to_string(lev) + " box " + std::to_string(i)
                     + " upper end not aligned to blocking factor in direction "
                     + std::to_string(d));
            }
        }
    }

    if (mesh.noChopDir() >= 0) {
        check_no_faces(ba, mesh.noChopDir(), lev);
        check_not_split(ba, mesh.noChopDir(), lev);
    }

    // Proper nesting relative to the coarser level.  All level lev grids
    // should be inside the grids on lev-1 with a buffer of
    // n_proper*bf_lev coarse cells, except at the physical boundary.
    {
        BoxArray const& cba = mesh.boxArray(lev-1);
        Box const& cdomain = mesh.Geom(lev-1).Domain();
        IntVect const np = mesh.bfLev(lev-1) * mesh.nProper();
        for (int i = 0; i < ba.size(); ++i) {
            Box b = amrex::coarsen(ba[i], rr);
            b.grow(np);
            b &= cdomain;
            if (!cba.contains(b, true)) {
                fail("level " + std::to_string(lev) + " box " + std::to_string(i)
                     + " not properly nested");
            }
        }
    }

    // Coverage: every tagged cell on lev-1 (plus its error buffer) must be
    // covered by the level lev grids, except cells outside the proper
    // nesting domain, which we skip for lev >= 2.
    if (lev == 1) {
        int const clev = lev-1;
        BoxArray const& cba = mesh.boxArray(clev);
        Box const& cdomain = mesh.Geom(clev).Domain();
        auto const& period = mesh.Geom(clev).periodicity();
        IntVect const nbuf = mesh.nErrorBufVect(clev);
        BoxArray const cfba = amrex::coarsen(ba, rr);
        iMultiFab mask = amrex::makeFineMask(cba, mesh.DistributionMap(clev),
                                             nbuf, ba, rr, period, 0, 1);
        Long ntagged = 0, nuncovered = 0;
        for (MFIter mfi(mask); mfi.isValid(); ++mfi) {
            auto const& m = mask.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k)
            {
                amrex::ignore_unused(k); // unused in 2D
                IntVect iv(AMREX_D_DECL(i,j,k));
                if (mesh.tagged(clev, iv)) {
                    ++ntagged;
                    Box nb(iv-nbuf, iv+nbuf);
                    bool ok = true;
                    amrex::LoopOnCpu(nb, [&] (int ii, int jj, int kk)
                    {
                        IntVect jv(AMREX_D_DECL(ii,jj,kk));
                        bool inside = true;
                        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                            if (!period.isPeriodic(d) &&
                                (jv[d] < cdomain.smallEnd(d) || jv[d] > cdomain.bigEnd(d))) {
                                inside = false;
                            }
                        }
                        if (inside && m(ii,jj,kk) == 0) { ok = false; }
                    });
                    if (!ok) { ++nuncovered; }
                }
            });
        }
        ParallelDescriptor::ReduceLongSum(ntagged);
        ParallelDescriptor::ReduceLongSum(nuncovered);
        amrex::Print() << "    tagged cells on level " << clev << ": " << ntagged
                       << ", uncovered: " << nuncovered << '\n';
        if (nuncovered > 0) { fail("tagged cells not covered by level 1 grids"); }
    }
}

// Linear function of position, constant in periodic directions.
struct LinearField
{
    Geometry geom;
    IndexType ityp;
    RealVect slope;

    [[nodiscard]] Real operator() (int i, int j, int k) const
    {
        amrex::ignore_unused(k); // unused in 2D
        IntVect const iv(AMREX_D_DECL(i,j,k));
        auto const dx = geom.CellSizeArray();
        auto const plo = geom.ProbLoArray();
        Real f = 1.0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            Real const off = ityp.nodeCentered(d) ? Real(0.0) : Real(0.5);
            f += slope[d] * (plo[d] + (iv[d]+off)*dx[d]);
        }
        return f;
    }

    void fillAll (MultiFab& mf) const
    {
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            auto const& a = mf.array(mfi);
            amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k)
            {
                a(i,j,k) = (*this)(i,j,k);
            });
        }
    }

    // Physical boundary functor for FillPatch: fill the cells outside the
    // domain (in non-periodic directions) with the exact values.
    void operator() (MultiFab& mf, int /*dcomp*/, int /*ncomp*/, IntVect const& /*nghost*/,
                     Real /*time*/, int /*bccomp*/) const
    {
        Box const pdomain = amrex::convert(geom.growPeriodicDomain(1024), ityp);
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            auto const& a = mf.array(mfi);
            amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k)
            {
                if (!pdomain.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                    a(i,j,k) = (*this)(i,j,k);
                }
            });
        }
    }
};

// FillPatchTwoLevels on the level lev grids, the way ERF calls it, with a
// linear field that must be reproduced exactly.  The face interpolater
// uses one-sided slopes next to non-periodic domain boundaries, so with
// interior_only the check skips two coarse cells next to those boundaries.
void test_fillpatch (TestMesh const& mesh, int lev, IndexType ityp,
                     Interpolater* mapper, std::string const& name, bool interior_only)
{
    Geometry const& cgeom = mesh.Geom(lev-1);
    Geometry const& fgeom = mesh.Geom(lev);
    IntVect const rr = mesh.refRatio(lev-1);
    BoxArray const cba = amrex::convert(mesh.boxArray(lev-1), ityp);
    BoxArray const fba = amrex::convert(mesh.boxArray(lev), ityp);

    RealVect slope;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        slope[d] = fgeom.isPeriodic(d) ? Real(0.0) : Real(0.7) + Real(0.3)*d;
    }
    LinearField cfield{cgeom, ityp, slope};
    LinearField ffield{fgeom, ityp, slope};

    // The coarse data must cover the coarsened fine ghost region plus one
    // cell for the interpolation stencil.  Level 0 covers the whole domain;
    // on finer coarse levels the proper nesting buffer limits the ghosts.
    IntVect const ngc(6);
    IntVect ngf(4);
    if (lev > 1) {
        IntVect const np = mesh.bfLev(lev-1) * mesh.nProper();
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            ngf[d] = std::clamp((np[d]-1)*rr[d], 0, 4);
        }
    }
    MultiFab cmf(cba, mesh.DistributionMap(lev-1), 1, ngc);
    MultiFab fmf(fba, mesh.DistributionMap(lev), 1, ngf);
    MultiFab dst(fba, mesh.DistributionMap(lev), 1, ngf);
    cfield.fillAll(cmf);
    ffield.fillAll(fmf);
    dst.setVal(Real(-1.e30));

    Vector<BCRec> bcs(1);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        int const bct = fgeom.isPeriodic(d) ? BCType::int_dir : BCType::foextrap;
        bcs[0].setLo(d, bct);
        bcs[0].setHi(d, bct);
    }

    Vector<MultiFab*> cmfv{&cmf};
    Vector<MultiFab*> fmfv{&fmf};
    Vector<Real> tv{Real(0.0)};
    amrex::FillPatchTwoLevels(dst, ngf, Real(0.0), cmfv, tv, fmfv, tv, 0, 0, 1,
                              cgeom, fgeom, cfield, 0, ffield, 0, rr, mapper, bcs, 0);

    Box const fdomain = amrex::convert(fgeom.Domain(), ityp);
    Box check_region = amrex::grow(fdomain, 1024);
    if (interior_only) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (!fgeom.isPeriodic(d)) {
                check_region.setSmall(d, fdomain.smallEnd(d) + 2*rr[d]);
                check_region.setBig  (d, fdomain.bigEnd  (d) - 2*rr[d]);
            }
        }
    }
    Real maxerr = 0;
    for (MFIter mfi(dst); mfi.isValid(); ++mfi) {
        auto const& a = dst.const_array(mfi);
        amrex::LoopOnCpu(mfi.fabbox() & check_region, [&] (int i, int j, int k)
        {
            maxerr = std::max(maxerr, std::abs(a(i,j,k) - ffield(i,j,k)));
        });
    }
    ParallelDescriptor::ReduceRealMax(maxerr);
    amrex::Print() << "    FillPatchTwoLevels " << name << " on level " << lev
                   << ": max error " << maxerr << '\n';
    if (maxerr > Real(1.e-10)) {
        fail("FillPatchTwoLevels " + name + " on level " + std::to_string(lev)
             + " is not exact for a linear field");
    }
}

// Coarsening a TagBoxArray whose grids are not aligned to the coarsening
// ratio produces overlapping valid regions.  Make sure the tags are still
// collated without duplicates.
void test_tag_overlap ()
{
    Box domain(IntVect(0), IntVect(AMREX_D_DECL(148,130,20)));
    BoxArray ba(domain);
    ba.maxSize(37);
    DistributionMapping dm(ba);
    IntVect const ratio(AMREX_D_DECL(8,8,8));
    AMREX_ALWAYS_ASSERT(!ba.coarsenable(ratio));

    for (auto const& pshift : {IntVect(0), IntVect(1)}) {
        Array<int,AMREX_SPACEDIM> is_per;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { is_per[d] = pshift[d]; }
        // ghost cells to exercise the buffer/ghost path as well
        TagBoxArray tags(ba, dm, IntVect(2));
        tags.setVal(ba, TagBox::SET);
        tags.buffer(IntVect(2));
        tags.coarsen(ratio);
        Box cdomain = amrex::coarsen(domain, ratio);
        Geometry cgeom(cdomain, RealBox(AMREX_D_DECL(0.,0.,0.), AMREX_D_DECL(1.,1.,1.)),
                       0, is_per);
        tags.mapPeriodicRemoveDuplicates(cgeom);
        Gpu::PinnedVector<IntVect> v;
        tags.collate(v);
        if (ParallelDescriptor::IOProcessor()) {
            Long const n = Long(v.size());
            std::sort(v.begin(), v.end());
            Long const nunique = std::unique(v.begin(), v.end()) - v.begin();
            // Buffered tags that fall outside a non-periodic boundary are
            // collated too, so only the duplicates are checked.
            amrex::Print() << "Tag overlap test (periodic " << pshift << "): "
                           << n << " tags, " << nunique << " unique, domain has "
                           << cdomain.numPts() << " cells\n";
            if (n != nunique) { fail("duplicate tags after coarsening"); }
            Long ninside = 0;
            for (auto const& iv : v) { if (cdomain.contains(iv)) { ++ninside; } }
            if (ninside != cdomain.numPts()) { fail("coarsened tags do not cover the domain"); }
        }

        // A single tag in a fine cell whose coarsened cell is also covered
        // by the valid region of another (untagged) fab must survive, with
        // and without ghost cells.  (Periodic images outside the domain
        // are collated too, so only the tags inside the domain count.)
        for (int ng = 0; ng <= 2; ng += 2) {
            TagBoxArray tags1(ba, dm, IntVect(ng));
            IntVect const iv1(AMREX_D_DECL(37,38,0));
            for (MFIter mfi(tags1); mfi.isValid(); ++mfi) {
                if (mfi.validbox().contains(iv1)) {
                    tags1.array(mfi)(AMREX_D_DECL(iv1[0],iv1[1],iv1[2])) = TagBox::SET;
                }
            }
            tags1.coarsen(ratio);
            tags1.mapPeriodicRemoveDuplicates(cgeom);
            Gpu::PinnedVector<IntVect> v1;
            tags1.collate(v1);
            if (ParallelDescriptor::IOProcessor()) {
                Long ninside = 0;
                bool right_cell = true;
                for (auto const& iv : v1) {
                    if (cdomain.contains(iv)) {
                        ++ninside;
                        if (iv != amrex::coarsen(iv1, ratio)) { right_cell = false; }
                    }
                }
                amrex::Print() << "Single tag test (periodic " << pshift << ", ngrow " << ng
                               << "): " << ninside << " tags inside the domain\n";
                if (ninside != 1 || !right_cell) {
                    fail("single tag lost or duplicated after coarsening");
                }
            }
        }
    }
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        // The mesh must be constructed first: it sets up the default
        // Geometry (including periodicity) from the inputs.
        TestMesh mesh;
        ParmParse pp("test");

        test_tag_overlap();

        {
            Vector<int> ebf;
            if (pp.queryarr("expected_bf_lev", ebf)) {
                IntVect const bf = mesh.bfLev(0);
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    if (bf[d] != ebf[d]) {
                        fail("bfLev(0) is " + std::to_string(bf[d]) + " in direction "
                             + std::to_string(d) + ", expected " + std::to_string(ebf[d]));
                    }
                }
            }
        }

        // Optionally chop the grids as if there were this many processes,
        // to exercise refine_grid_layout deterministically.
        int chop_target = 0;
        pp.query("chop_target", chop_target);

        BoxArray ba0 = mesh.MakeBaseGrids();
        if (chop_target > 0) { mesh.ChopGrids(0, ba0, chop_target); }
        mesh.SetBoxArray(0, ba0);
        mesh.SetDistributionMap(0, DistributionMapping(ba0));
        mesh.SetFinestLevel(0);
        amrex::Print() << "Level 0: " << ba0.size() << " grids\n";
        check_level0(mesh);

        Vector<BoxArray> new_grids;
        int new_finest = 0;
        // Add one level at a time until all levels exist, regridding from
        // level 0 each time.
        for (int iter = 0; iter < mesh.maxLevel()+1; ++iter) {
            int const lbase = 0;
            auto const t0 = amrex::second();
            mesh.MakeNewGrids(lbase, 0.0, new_finest, new_grids);
            amrex::Print() << "MakeNewGrids from level 0 took " << amrex::second()-t0 << " s\n";
            for (int lev = lbase+1; lev <= new_finest; ++lev) {
                if (chop_target > 0) {
                    Long const nbefore = new_grids[lev].size();
                    mesh.ChopGrids(lev, new_grids[lev], chop_target);
                    amrex::Print() << "    ChopGrids on level " << lev << ": " << nbefore
                                   << " -> " << new_grids[lev].size() << " grids\n";
                }
                mesh.SetBoxArray(lev, new_grids[lev]);
                mesh.SetDistributionMap(lev, DistributionMapping(new_grids[lev]));
            }
            mesh.SetFinestLevel(new_finest);
        }
        // Regrid from the finest base level too.
        if (mesh.finestLevel() >= 2) {
            int const lbase = mesh.finestLevel()-1;
            mesh.MakeNewGrids(lbase, 0.0, new_finest, new_grids);
            for (int lev = lbase+1; lev <= new_finest; ++lev) {
                if (chop_target > 0) {
                    Long const nbefore = new_grids[lev].size();
                    mesh.ChopGrids(lev, new_grids[lev], chop_target);
                    amrex::Print() << "    ChopGrids on level " << lev << ": " << nbefore
                                   << " -> " << new_grids[lev].size() << " grids\n";
                }
                mesh.SetBoxArray(lev, new_grids[lev]);
                mesh.SetDistributionMap(lev, DistributionMapping(new_grids[lev]));
            }
            mesh.SetFinestLevel(new_finest);
        }

        int expected_finest = mesh.maxLevel();
        pp.query("expected_finest", expected_finest);
        if (mesh.finestLevel() != expected_finest) {
            fail("finest level is " + std::to_string(mesh.finestLevel())
                 + ", expected " + std::to_string(expected_finest));
        }

        for (int lev = 1; lev <= mesh.finestLevel(); ++lev) {
            check_level(mesh, lev);
            test_fillpatch(mesh, lev, IndexType::TheCellType(), &cell_cons_interp, "cell", false);
            test_fillpatch(mesh, lev, IndexType(IntVect::TheDimensionVector(0)),
                           &face_cons_linear_interp, "x-face", true);
        }

        std::string dump;
        pp.query("dump_grids", dump);
        if (!dump.empty() && ParallelDescriptor::IOProcessor()) {
            std::ofstream ofs(dump);
            for (int lev = 0; lev <= mesh.finestLevel(); ++lev) {
                ofs << "Level " << lev << '\n';
                for (int i = 0; i < mesh.boxArray(lev).size(); ++i) {
                    ofs << mesh.boxArray(lev)[i] << '\n';
                }
            }
        }

        ParallelDescriptor::ReduceIntSum(nfail);
        if (nfail == 0) {
            amrex::Print() << "PASS\n";
        } else {
            amrex::Abort("Grid creation test failed");
        }
    }
    amrex::Finalize();
}
