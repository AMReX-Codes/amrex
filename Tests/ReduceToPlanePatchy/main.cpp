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
        DistributionMapping dm(ba);
        MultiFab mf(ba, dm, 1, 0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
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
        ParallelDescriptor::ReduceRealSum(ref_plane.dataPtr(),
                                          static_cast<int>(ref_plane.box().numPts()));

        auto [plane_patch, plane_unique] = ReduceToPlaneMF2Patchy<ReduceOpSum>(dir, domain, mf,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) -> Real
            {
                return ma[box_no](i,j,k);
            });

        // Sanity: one projected FAB per input FAB.
        AMREX_ALWAYS_ASSERT(plane_patch.size() == mf.size());

        Gpu::streamSynchronize();

        // Compare unique sparse result to reference plane on overlapping cells.
        Real max_err = 0.0;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion()) reduction(max:max_err)
#endif
        for (MFIter mfi(plane_unique); mfi.isValid(); ++mfi) {
            auto const& pa = plane_unique.const_array(mfi);
            Box const& bx = mfi.validbox();
            Real local_max = 0.0;
            AMREX_LOOP_3D(bx, i, j, k,
            {
                local_max = std::max(local_max, std::abs(pa(i,j,k) - ref_plane(IntVect(AMREX_D_DECL(i, j, k)))));
            });
            max_err = std::max(max_err, local_max);
        }

        AMREX_ALWAYS_ASSERT(max_err == 0.0);
    }

    amrex::Finalize();
    return 0;
}
