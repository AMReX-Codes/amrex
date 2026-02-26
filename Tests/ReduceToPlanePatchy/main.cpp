#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <algorithm>
#include <cmath>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        constexpr int dir = 0;
        const Box domain(IntVect(0), IntVect(AMREX_D_DECL(7, 7, 7)));

        // Patchy BoxArray: leaves a gap in x so this is not a rectangular domain coverage.
        BoxList bl;
        bl.push_back(Box(IntVect(AMREX_D_DECL(0, 0, 0)), IntVect(AMREX_D_DECL(1, 7, 7))));
        bl.push_back(Box(IntVect(AMREX_D_DECL(4, 0, 0)), IntVect(AMREX_D_DECL(5, 7, 7))));
        bl.push_back(Box(IntVect(AMREX_D_DECL(6, 0, 0)), IntVect(AMREX_D_DECL(7, 7, 7))));

        BoxArray ba(std::move(bl));
        ba.maxSize(4);
        DistributionMapping dm(ba);
        MultiFab mf(ba, dm, 1, 0);

        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            auto const& a = mf.array(mfi);
            Box const& bx = mfi.validbox();
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                a(i,j,k) = static_cast<Real>(1 + i + 10*j + 100*k);
            });
        }

        auto const& ma = mf.const_arrays();

        // Reference: BaseFab plane reduction (works on patchy layouts).
        auto ref_plane = ReduceToPlane<ReduceOpSum,Real>(dir, domain, mf,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) -> Real
            {
                return ma[box_no](i,j,k);
            });
        if (ParallelDescriptor::NProcs() > 1) {
            auto npts = static_cast<int>(ref_plane.box().numPts());
#ifdef AMREX_USE_GPU
            auto* dp = ref_plane.dataPtr();
            Gpu::PinnedVector<Real> hv(npts);
            auto* hp = hv.data();
            Gpu::copyAsync(Gpu::deviceToHost, dp, dp+npts, hp);
            Gpu::streamSynchronize();
#else
            auto* hp = ref_plane.dataPtr();
#endif
            ParallelDescriptor::ReduceRealSum(hp, npts);
#ifdef AMREX_USE_GPU
            Gpu::copyAsync(Gpu::hostToDevice, hp, hp+npts, dp);
            Gpu::streamSynchronize();
#endif
        }

        auto [plane_patch, plane_unique] = ReduceToPlaneMF2Patchy<ReduceOpSum>(dir, domain, mf,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) -> Real
            {
                return ma[box_no](i,j,k);
            });

        // Sanity: one projected FAB per input FAB.
        AMREX_ALWAYS_ASSERT(plane_patch.size() == mf.size());

        Gpu::streamSynchronize();

        // Compare unique sparse result to reference plane on overlapping cells.
        auto const& res = plane_unique.const_arrays();
        auto const& ref = ref_plane.const_array();
        Real max_err = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{},
                                 plane_unique,
                                 [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
                                     -> GpuTuple<Real>
                                 {
                                     return {std::abs(res[b](i,j,k) - ref(i,j,k))};
                                 });
        ParallelDescriptor::ReduceRealMax(max_err);
        AMREX_ALWAYS_ASSERT(amrex::almostEqual(max_err, Real(0.0)));
    }

    amrex::Finalize();
    return 0;
}
