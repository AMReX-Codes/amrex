#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_OpenBC.H>

#include <cmath>

using namespace amrex;

namespace {

void solve_fully_refined_hierarchy (int nlevels)
{
    constexpr int n_cell = 8;

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

        grids[lev].define(domain);
        grids[lev].maxSize(16);
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
        solve_fully_refined_hierarchy(2);
        solve_fully_refined_hierarchy(3);
        solve_fully_refined_hierarchy(4);
        solve_fully_refined_hierarchy(5);
    }
    amrex::Finalize();
}
