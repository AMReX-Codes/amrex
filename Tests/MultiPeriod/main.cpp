#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_FabArray.H>
#include <AMReX_MFParallelFor.H>
#include <AMReX_ParReduce.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        // Domain size: 2 x 128 x 4
        Box box(IntVect(0), IntVect(AMREX_D_DECL(1, 127, 3)));
        Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,1)};
        Geometry geom(box, RealBox(AMREX_D_DECL(Real(0),Real(0),Real(0)),
                                   AMREX_D_DECL(Real(1),Real(1),Real(1))),
                      CoordSys::cartesian, is_periodic);
        BoxArray ba(box);
        ba.maxSize(32);
        ba.convert(IntVect(AMREX_D_DECL(1,0,0))); // nodal in x-direction
        DistributionMapping dm(ba);

        FabArray<BaseFab<Long>> mf1(ba,dm,1,IntVect(4));
        FabArray<BaseFab<Long>> mf2(ba,dm,1,IntVect(5));

        mf1.setVal(-1);
        mf2.setVal(-2);

        auto const& len = geom.Domain().length3d();
        auto expected = [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                while (i <  0     ) { i += len[0]; }
                while (i >= len[0]) { i -= len[0]; }
                while (j <  0     ) { j += len[1]; }
                while (j >= len[1]) { j -= len[1]; }
                while (k <  0     ) { k += len[2]; }
                while (k >= len[2]) { k -= len[2]; }
                return Long(i) + Long(j)*Long(len[0]) + Long(k)*Long(len[0])*Long(len[1]);
            };

        auto const& ma1 = mf1.arrays();
        auto const& ma2 = mf2.arrays();

        // Initialize valid region
        ParallelFor(mf1, IntVect(0), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            ma1[b](i,j,k) = expected(i,j,k);
        });

        mf1.FillBoundary(geom.periodicity());
        mf2.ParallelCopy(mf1, 0, 0, 1, IntVect(0), mf2.nGrowVect(), geom.periodicity());

        auto r1 = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Long>{}, mf1, mf1.nGrowVect(),
                            [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) -> GpuTuple<Long>
                            {
                                return { Long(expected(i,j,k) != ma1[b](i,j,k)) };
                            });
        auto r2 = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Long>{}, mf2, mf2.nGrowVect(),
                            [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) -> GpuTuple<Long>
                            {
                                return { Long(expected(i,j,k) != ma2[b](i,j,k)) };
                            });

        AMREX_ALWAYS_ASSERT(r1 == 0);
        AMREX_ALWAYS_ASSERT(r2 == 0);

        if (r1 == 0 && r2 == 0) {
            amrex::Print() << "SUCCESS\n";
        }
    }
    amrex::Finalize();
}
