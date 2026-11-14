#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_OpenBC.H>

#include <cmath>

using namespace amrex;

namespace {

enum struct Refinement { full, corner, face };

void solve_hierarchy (int nlevels, Refinement refinement)
{
    // full: every level covers its whole domain.
    // corner: every fine level is a single 32^3 grid touching the low corner
    // of the domain, so support grids are needed on many levels and along
    // the exterior edges and corner.
    // face: every fine level is an 8x16x16 patch of 8^3 grids in the middle
    // of the high-x face. Each patch is nested in its parent with a margin of
    // only 4 coarse cells along the face, so support grids must not be wider
    // than the grids they extend.
    int const n_cell = (refinement == Refinement::full) ? 8 : 32;

    Vector<Geometry> geom(nlevels);
    Vector<BoxArray> grids(nlevels);
    Vector<DistributionMapping> dmap(nlevels);
    Vector<MultiFab> phi(nlevels);
    Vector<MultiFab> rhs(nlevels);

    RealBox const real_box(0._rt, 0._rt, 0._rt, 1._rt, 1._rt, 1._rt);
    Array<int, AMREX_SPACEDIM> const periodic{0, 0, 0};

    for (int lev = 0; lev < nlevels; ++lev) {
        int const nc = n_cell * (1 << lev);
        Box const domain(IntVect(0), IntVect(nc-1));
        geom[lev].define(domain, real_box, CoordSys::cartesian, periodic);

        if (refinement == Refinement::corner && lev > 0) {
            grids[lev].define(Box(IntVect(0), IntVect(31)));
        } else if (refinement == Refinement::face && lev > 0) {
            grids[lev].define(Box(IntVect(nc-8, nc/2-8, nc/2-8),
                                  IntVect(nc-1, nc/2+7, nc/2+7)));
        } else {
            grids[lev].define(domain);
        }
        grids[lev].maxSize((refinement == Refinement::face) ? 8 : 16);
        dmap[lev].define(grids[lev]);

        phi[lev].define(grids[lev], dmap[lev], 1, 1);
        rhs[lev].define(grids[lev], dmap[lev], 1, 2);
        phi[lev].setVal(0._rt);
        rhs[lev].setVal(0._rt);

        auto const dx = geom[lev].CellSizeArray();
        auto const problo = geom[lev].ProbLoArray();
        for (MFIter mfi(rhs[lev]); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.validbox();
            auto const& a = rhs[lev].array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real const x = problo[0] + (i+0.5_rt)*dx[0] - 0.5_rt;
                Real const y = problo[1] + (j+0.5_rt)*dx[1] - 0.5_rt;
                Real const z = problo[2] + (k+0.5_rt)*dx[2] - 0.5_rt;
                a(i,j,k) = std::exp(-64._rt*(x*x+y*y+z*z));
            });
        }
    }

    LPInfo info;
    info.setDeterministic(true);
    OpenBCSolver solver(geom, grids, dmap, info);
    Real const error = solver.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs),
                                    1.e-10_rt, 0._rt);
    AMREX_ALWAYS_ASSERT(std::isfinite(error));

    for (auto const& mf : phi) {
        AMREX_ALWAYS_ASSERT(!mf.contains_nan());
        AMREX_ALWAYS_ASSERT(!mf.contains_inf());
    }
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        for (int nlevels = 2; nlevels <= 5; ++nlevels) {
            solve_hierarchy(nlevels, Refinement::full);
        }
        for (int nlevels = 2; nlevels <= 7; ++nlevels) {
            solve_hierarchy(nlevels, Refinement::corner);
        }
        for (int nlevels = 2; nlevels <= 7; ++nlevels) {
            solve_hierarchy(nlevels, Refinement::face);
        }
    }
    amrex::Finalize();
}
