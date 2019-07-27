
#include <AMReX.H>
#include <AMReX_Utility.H>
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
    int ncell = 512;
    int max_grid_size = 128;
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
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.validbox();
        auto const& fab = mf.array(mfi);
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fab(i,j,k) = amrex::Random();
        });
    }

    {
        BL_PROFILE("MultiFab::sum");
        amrex::Print().SetPrecision(17) << "sum: " << mf.sum() << "\n";
    }

    {
        BL_PROFILE("MultiFab::min");
        amrex::Print().SetPrecision(17) << "min: " << mf.min(0) << "\n";
    }

    {
        BL_PROFILE("MultiFab::max");
        amrex::Print().SetPrecision(17) << "max: " << mf.max(0) << "\n";
    }

    {
        BL_PROFILE("ParallelReduce3");
        Vector<Real> hv{               0.0, // initial value for sum
                std::numeric_limits<Real>::max(),             // min
                std::numeric_limits<Real>::lowest()};         // max
        Gpu::DeviceVector<Real> dv(hv);
        Real* dp = dv.dataPtr();

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            auto const& fab = mf.array(mfi);
            typedef GpuTuple<Real,Real,Real> Real3;
            amrex::ParallelForReduce
                (bx, // Box
                 Real3{hv[0],hv[1],hv[2]}, // initial values
            // First lambda works on each cell
            [=] AMREX_GPU_DEVICE (int i, int j, int k, Real3& r)  noexcept
            {
                Real x = fab(i,j,k);
                Real& r0 = amrex::get<0>(r);
                Real& r1 = amrex::get<1>(r);
                Real& r2 = amrex::get<2>(r);
                r0 += x;
                r1 = amrex::min(r1,x);
                r2 = amrex::max(r2,x);
            },
            // Second lambda does reduce and stores the result in global memory
            [=] AMREX_GPU_DEVICE (Real3 const& r) noexcept
            {
                amrex::ReduceSum(dp  , amrex::get<0>(r));
                amrex::ReduceMin(dp+1, amrex::get<1>(r));
                amrex::ReduceMax(dp+2, amrex::get<2>(r));
            });
        }

        Gpu::dtoh_memcpy(hv.data(), dp, hv.size()*sizeof(Real));
        // MPI reduce
        ParallelDescriptor::ReduceRealSum(hv[0]);
        ParallelDescriptor::ReduceRealMin(hv[1]);
        ParallelDescriptor::ReduceRealMax(hv[2]);
        amrex::Print().SetPrecision(17) << "sum: " << hv[0] << "\n"
                                        << "min: " << hv[1] << "\n"
                                        << "max: " << hv[2] << "\n";
    }
}
