#include <AMReX.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Reduce.H>

using namespace amrex;

// Two interpolaters with different CoarseBox() share the same FillPatch
// metadata cache entry. The second fill must still see coarse patches
// large enough for its stencil.

namespace {

void fill_linear (MultiFab& mf, Geometry const& geom, IntVect const& ng)
{
    auto const& dx = geom.CellSizeArray();
    auto const& ma = mf.arrays();
    ParallelFor(mf, ng, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
    {
        amrex::ignore_unused(j,k);
        ma[b](i,j,k) = AMREX_D_TERM(  (Real(i)+Real(0.5))*dx[0],
                                    + Real(2.)*(Real(j)+Real(0.5))*dx[1],
                                    + Real(3.)*(Real(k)+Real(0.5))*dx[2]);
    });
    Gpu::streamSynchronize();
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        int const ncell = 32;    // coarse cells per direction
        int const nghost = 3;    // odd: grown patch boxes are misaligned
        IntVect const ratio(2);

        Box const cdom(IntVect(0), IntVect(ncell-1));
        Box const fdom = amrex::refine(cdom, ratio);
        RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
        Array<int,AMREX_SPACEDIM> isper{AMREX_D_DECL(0,0,0)};
        Geometry const cgeom(cdom, rb, CoordSys::cartesian, isper);
        Geometry const fgeom(fdom, rb, CoordSys::cartesian, isper);

        BoxArray cba(cdom); cba.maxSize(16);
        DistributionMapping cdm(cba);
        BoxArray fba(Box(IntVect(ncell/2), IntVect(ncell/2+ncell-1))); fba.maxSize(16);
        DistributionMapping fdm(fba);

        Vector<BCRec> bcs(1);
        for (OrientationIter oit; oit; ++oit) { bcs[0].set(oit(), BCType::foextrap); }
        PhysBCFunct<GpuBndryFuncFab<FabFillNoOp>> cbf(cgeom, bcs, GpuBndryFuncFab<FabFillNoOp>{});
        PhysBCFunct<GpuBndryFuncFab<FabFillNoOp>> fbf(fgeom, bcs, GpuBndryFuncFab<FabFillNoOp>{});

        MultiFab cmf(cba, cdm, 1, 1);
        MultiFab fmf(fba, fdm, 1, nghost);
        fill_linear(cmf, cgeom, IntVect(1));
        fill_linear(fmf, fgeom, IntVect(0));

        // First fill builds the cache entry with the smaller bilinear coarse boxes.
        {
            MultiFab tmp(fba, fdm, 1, nghost);
            tmp.setVal(0.);
            InterpBase* mapper = &mf_cell_bilinear_interp;
            FillPatchTwoLevels(tmp, IntVect(nghost), Real(0.), {&cmf}, {Real(0.)},
                               {&fmf}, {Real(0.)}, 0, 0, 1, cgeom, fgeom,
                               cbf, 0, fbf, 0, ratio, mapper, bcs, 0);
        }

        // Second fill needs coarse boxes grown by one on every side.
        {
            InterpBase* mapper = &mf_lincc_interp;
            FillPatchTwoLevels(fmf, IntVect(nghost), Real(0.), {&cmf}, {Real(0.)},
                               {&fmf}, {Real(0.)}, 0, 0, 1, cgeom, fgeom,
                               cbf, 0, fbf, 0, ratio, mapper, bcs, 0);
        }

        // Linear data is reproduced exactly inside the domain.
        auto const& dx = fgeom.CellSizeArray();
        auto const& ma = fmf.const_arrays();
        Real err = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, fmf, IntVect(nghost),
            [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept -> GpuTuple<Real>
            {
                amrex::ignore_unused(j,k);
                if (!fdom.contains(i,j,k)) { return Real(0.); }
                Real ex = AMREX_D_TERM(  (Real(i)+Real(0.5))*dx[0],
                                       + Real(2.)*(Real(j)+Real(0.5))*dx[1],
                                       + Real(3.)*(Real(k)+Real(0.5))*dx[2]);
                return std::abs(ma[b](i,j,k) - ex);
            });
        ParallelDescriptor::ReduceRealMax(err);
        amrex::Print() << "FillPatch cache test: max error " << err << "\n";
        Real const tol = Real(100.)*std::numeric_limits<Real>::epsilon();
        if (!(err < tol)) {
            amrex::Abort("FillPatch cache test failed");
        }
    }
    amrex::Finalize();
}
