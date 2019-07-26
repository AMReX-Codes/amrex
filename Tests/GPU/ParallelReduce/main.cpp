
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

void main_main();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    main_main();
    amrex::Finalize();
}


void main_main ()
{
    int ncell = 256;
    int max_grid_size = 64;
    {
        ParmParse pp;
        pp.query("ncell", ncell);
        pp.query("max_grid_size", max_grid_size);
    }

    BoxArray ba;
    {
        Box domain_box(IntVect(0), IntVect(ncell-1));
        ba.define(domain_box);
        ba.maxSize(max_grid_size);
    }

    MultiFab mf(ba,DistributionMapping{ba},1,0);
    mf.setVal(1.0);

    {
        BL_PROFILE("MultiFab::sum");
        amrex::Print() << mf.sum() << "\n";
    }

    {
        BL_PROFILE("shfl_down");
        Gpu::DeviceScalar<Real> r(0.0);
        Real* dp = r.dataPtr();
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto const& arr = mf.array(mfi);

            const auto ec = amrex::Gpu::ExecutionConfig(bx.numPts()/32);
            amrex::launch_global<<<ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream()>>>(
            [=] AMREX_GPU_DEVICE ()
            {
                const BaseFab<Real> fab(arr,bx.ixType());
                Real tsum = 0.0;
                for (auto const tbx : Gpu::Range(bx)) {
                    tsum += fab.sum(tbx, 0);
                }
                tsum = Gpu::warpReduceSum<32,Real>(tsum);
                // if (threadIdx.x % 32 == 0) {
                if ((threadIdx.x & 0x1f) == 0) {
                    Gpu::Atomic::Add(dp, tsum);
                }
            });

        }
        amrex::Print() << r.dataValue() << "\n";
    }

    {
        BL_PROFILE("MultiFab::sum_22");
        amrex::Print() << mf.sum() << "\n";
    }

    {
        BL_PROFILE("shfl_down_shared");
        Gpu::DeviceScalar<Real> r(0.0);
        Real* dp = r.dataPtr();
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto const& arr = mf.array(mfi);

            const auto ec = amrex::Gpu::ExecutionConfig(bx.numPts()/32);
            amrex::launch_global<<<ec.numBlocks, ec.numThreads, 32*sizeof(Real), Gpu::gpuStream()>>>(
            [=] AMREX_GPU_DEVICE ()
            {
                const BaseFab<Real> fab(arr,bx.ixType());
                Real tsum = 0.0;
                for (auto const tbx : Gpu::Range(bx)) {
                    tsum += fab.sum(tbx, 0);
                }
                tsum = Gpu::blockReduceSum<32,Real>(tsum);
                // if (threadIdx.x % 32 == 0) {
                if (threadIdx.x == 0) {
                    Gpu::Atomic::Add(dp, tsum);
                }
            });

        }
        amrex::Print() << r.dataValue() << "\n";
    }

}
