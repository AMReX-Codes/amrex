#include <AMReX.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Reduce.H>

#include <cmath>
#include <limits>

using namespace amrex;

// Two fills of the same fine MultiFab with different interpolaters or
// ratios look up the same FillPatch metadata cache key. The second fill
// must still see coarse patches large enough for its stencil.

namespace {

void fill_linear (MultiFab& mf, Geometry const& geom, IntVect const& ng)
{
    auto const& dx = geom.CellSizeArray();
    Real const off = mf.ixType().cellCentered() ? Real(0.5) : Real(0.);
    auto const& ma = mf.arrays();
    ParallelFor(mf, ng, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
    {
        amrex::ignore_unused(j,k);
        ma[b](i,j,k) = AMREX_D_TERM(  (Real(i)+off)*dx[0],
                                    + Real(2.)*(Real(j)+off)*dx[1],
                                    + Real(3.)*(Real(k)+off)*dx[2]);
    });
    Gpu::streamSynchronize();
}

// Fill with `first` from a coarse level at ratio `first_ratio`, then with
// `second` at `ratio` on the same fine MultiFab, and return the max error
// of the second fill inside the domain.
Real run (IntVect const& ratio, int nghost, bool periodic, InterpBase* first, InterpBase* second,
          IndexType typ = IndexType::TheCellType(), IntVect first_ratio = IntVect(0))
{
    if (first_ratio == IntVect(0)) { first_ratio = ratio; }

    int const ncell = 32;    // coarse cells per direction

    Box const cdom(IntVect(0), IntVect(ncell-1));
    Box const fdom = amrex::refine(cdom, ratio);
    RealBox rb(AMREX_D_DECL(0.,0.,0.), AMREX_D_DECL(1.,1.,1.));
    Array<int,AMREX_SPACEDIM> isper{AMREX_D_DECL(periodic,periodic,periodic)};
    Geometry const cgeom(cdom, rb, CoordSys::cartesian, isper);
    Geometry const fgeom(fdom, rb, CoordSys::cartesian, isper);

    Box const cdom1 = amrex::coarsen(fdom, first_ratio);
    Geometry const cgeom1(cdom1, rb, CoordSys::cartesian, isper);

    BoxArray cba(cdom); cba.maxSize(16);
    DistributionMapping cdm(cba);
    BoxArray cba1(cdom1); cba1.maxSize(16);
    DistributionMapping cdm1(cba1);
    // Fine level: checkerboard of 4^d blocks well inside the domain, so
    // periodic images are never used and FPinfo sees more than 8 boxes
    // after simplification in 2D and 3D.
    BoxArray fba(Box(IntVect(ncell/2), IntVect(ncell/2+ncell-1))); fba.maxSize(4);
    BoxList fbl;
    for (int i = 0; i < static_cast<int>(fba.size()); ++i) {
        Box const b = fba[i];
        if ((b.smallEnd()/4).sum() % 2 == 0) { fbl.push_back(b); }
    }
    fba = BoxArray(std::move(fbl));
    DistributionMapping fdm(fba);
    cba.convert(typ);
    cba1.convert(typ);
    fba.convert(typ);

    Vector<BCRec> bcs(1);
    for (OrientationIter oit; oit; ++oit) { bcs[0].set(oit(), BCType::foextrap); }
    PhysBCFunct<GpuBndryFuncFab<FabFillNoOp>> cbf(cgeom, bcs, GpuBndryFuncFab<FabFillNoOp>{});
    PhysBCFunct<GpuBndryFuncFab<FabFillNoOp>> cbf1(cgeom1, bcs, GpuBndryFuncFab<FabFillNoOp>{});
    PhysBCFunct<GpuBndryFuncFab<FabFillNoOp>> fbf(fgeom, bcs, GpuBndryFuncFab<FabFillNoOp>{});

    MultiFab cmf(cba, cdm, 1, 1);
    MultiFab cmf1(cba1, cdm1, 1, 1);
    MultiFab fmf(fba, fdm, 1, nghost);
    fill_linear(cmf, cgeom, IntVect(1));
    fill_linear(cmf1, cgeom1, IntVect(1));
    fill_linear(fmf, fgeom, IntVect(0));

    {
        MultiFab tmp(fba, fdm, 1, nghost);
        tmp.setVal(0.);
        FillPatchTwoLevels(tmp, IntVect(nghost), Real(0.), {&cmf1}, {Real(0.)},
                           {&fmf}, {Real(0.)}, 0, 0, 1, cgeom1, fgeom,
                           cbf1, 0, fbf, 0, first_ratio, first, bcs, 0);
    }
    FillPatchTwoLevels(fmf, IntVect(nghost), Real(0.), {&cmf}, {Real(0.)},
                       {&fmf}, {Real(0.)}, 0, 0, 1, cgeom, fgeom,
                       cbf, 0, fbf, 0, ratio, second, bcs, 0);

    // Linear data is reproduced exactly inside the domain.
    Box const fdom_t = amrex::convert(fdom, typ);
    auto const& dx = fgeom.CellSizeArray();
    Real const off = typ.cellCentered() ? Real(0.5) : Real(0.);
    auto const& ma = fmf.const_arrays();
    Real err = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, fmf, IntVect(nghost),
        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept -> GpuTuple<Real>
        {
            amrex::ignore_unused(j,k);
            if (!fdom_t.contains(i,j,k)) { return Real(0.); }
            Real ex = AMREX_D_TERM(  (Real(i)+off)*dx[0],
                                   + Real(2.)*(Real(j)+off)*dx[1],
                                   + Real(3.)*(Real(k)+off)*dx[2]);
            return std::abs(ma[b](i,j,k) - ex);
        });
    ParallelDescriptor::ReduceRealMax(err);
    return err;
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        Real const tol = Real(100.)*std::numeric_limits<Real>::epsilon();
        // Odd ghost width: grown patch boxes are misaligned.
        Real err1 = run(IntVect(2), 3, false, &mf_cell_bilinear_interp, &mf_lincc_interp);
        // Periodic domain grown by the ghost width, ratio 4: bilinear looks
        // like piecewise constant on the domain box and its neighbors.
        Real err2 = run(IntVect(4), 1, true, &mf_pc_interp, &mf_cell_bilinear_interp);
        // Nodal data, one-node ghost slabs.
        Real err3 = run(IntVect(2), 1, false, &pc_interp, &node_bilinear_interp,
                        IndexType::TheNodeType());
        // Same interpolater, but the first fill comes from a ratio-4 level.
        Real err4 = run(IntVect(2), 1, false, &mf_cell_bilinear_interp, &mf_cell_bilinear_interp,
                        IndexType::TheCellType(), IntVect(4));
        amrex::Print() << "FillPatch cache test: max errors " << err1 << " " << err2
                       << " " << err3 << " " << err4 << "\n";
        if (!(err1 < tol && err2 < tol && err3 < tol && err4 < tol)) {
            amrex::Abort("FillPatch cache test failed");
        }
    }
    amrex::Finalize();
}
