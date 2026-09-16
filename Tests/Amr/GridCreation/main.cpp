// Exercise AmrMesh grid creation and check the properties of the resulting
// grids: max_grid_size, blocking factor alignment, coarsenability by
// ref_ratio, coverage of tagged cells and proper nesting.

#include <AMReX.H>
#include <AMReX_AmrMesh.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_TagBox.H>
#include <AMReX_iMultiFab.H>

#include <algorithm>
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
        auto const& geom = Geom(lev);
        auto const dx = geom.CellSizeArray();
        auto const plo = geom.ProbLoArray();
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
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                int rr = refRatio(lev-1)[d];
                g[d] = std::max(1, blockingFactor(lev)[d]/rr) * rr;
            }
        }
        return g;
    }

    [[nodiscard]] IntVect bfLev (int lev) const
    {
        IntVect b(1);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            b[d] = std::max(1, blockingFactor(lev+1)[d]/refRatio(lev)[d]);
        }
        return b;
    }

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

void check_level (TestMesh const& mesh, int lev)
{
    BoxArray const& ba = mesh.boxArray(lev);
    Box const& domain = mesh.Geom(lev).Domain();
    IntVect const rr = mesh.refRatio(lev-1);
    IntVect const unit = mesh.gridUnit(lev);
    IntVect const mgs = mesh.maxGridSize(lev);

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
    }
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        test_tag_overlap();

        TestMesh mesh;
        BoxArray ba0 = mesh.MakeBaseGrids();
        mesh.SetBoxArray(0, ba0);
        mesh.SetDistributionMap(0, DistributionMapping(ba0));
        mesh.SetFinestLevel(0);
        amrex::Print() << "Level 0: " << ba0.size() << " grids\n";

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
                mesh.SetBoxArray(lev, new_grids[lev]);
                mesh.SetDistributionMap(lev, DistributionMapping(new_grids[lev]));
            }
            mesh.SetFinestLevel(new_finest);
        }

        int expected_finest = mesh.maxLevel();
        ParmParse pp("test");
        pp.query("expected_finest", expected_finest);
        if (mesh.finestLevel() != expected_finest) {
            fail("finest level is " + std::to_string(mesh.finestLevel())
                 + ", expected " + std::to_string(expected_finest));
        }

        for (int lev = 1; lev <= mesh.finestLevel(); ++lev) {
            check_level(mesh, lev);
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
