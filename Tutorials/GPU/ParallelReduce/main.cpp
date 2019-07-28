
#include <AMReX.H>
#include <AMReX_Utility.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

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
    iMultiFab imf(ba,mf.DistributionMap(),1,0);
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.validbox();
        auto const& fab = mf.array(mfi);
        auto const& ifab = imf.array(mfi);
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fab(i,j,k) = amrex::Random();
            ifab(i,j,k) = (amrex::Random() > 0.5) ? 1 : 0;
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
        BL_PROFILE("iMultiFab::sum");
        amrex::Print() << "isum: " << imf.sum(0) << "\n";
    }

    {
        BL_PROFILE("FabReduceTuple");

        Vector<Real> hv{               0.0, // initial value for sum
                std::numeric_limits<Real>::max(),             // min
                std::numeric_limits<Real>::lowest()};         // max
        Gpu::DeviceVector<Real> dv(hv);
        Real* dp = dv.dataPtr();

        int isum = 0; // initialize to 0
        Gpu::DeviceScalar<int> ds(isum); 
        int * ip = ds.dataPtr();

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            auto const& fab = mf.array(mfi);
            auto const& ifab = imf.array(mfi);
            typedef GpuTuple<Real,Real,Real,int> ReduceTuple;
            amrex::FabReduce
                (bx, // Box
                 ReduceTuple{hv[0],hv[1],hv[2], isum}, // initial values
            // First lambda works on each cell
            [=] AMREX_GPU_DEVICE (int i, int j, int k, ReduceTuple* r)  noexcept
            {
                Real x = fab(i,j,k);
                Real& r0 = amrex::get<0>(*r);
                Real& r1 = amrex::get<1>(*r);
                Real& r2 = amrex::get<2>(*r);
                r0 += x;
                r1 = amrex::min(r1,x);
                r2 = amrex::max(r2,x);
                amrex::get<3>(*r) += ifab(i,j,k);
            },
            // Second lambda does reduce and stores the result in global memory
            // Reduce must be done with Gpu::Reduce[Sum|Min|Max].
            [=] AMREX_GPU_DEVICE (ReduceTuple const& r) noexcept
            {
                Gpu::ReduceSum(dp  , amrex::get<0>(r));
                Gpu::ReduceMin(dp+1, amrex::get<1>(r));
                Gpu::ReduceMax(dp+2, amrex::get<2>(r));
                Gpu::ReduceSum(ip  , amrex::get<3>(r));
            });
        }

        Gpu::dtoh_memcpy(hv.data(), dp, hv.size()*sizeof(Real));
        isum = ds.dataValue();
        // MPI reduce
        ParallelDescriptor::ReduceRealSum(hv[0]);
        ParallelDescriptor::ReduceRealMin(hv[1]);
        ParallelDescriptor::ReduceRealMax(hv[2]);
        ParallelDescriptor::ReduceIntSum(isum);
        amrex::Print().SetPrecision(17) << "sum: " << hv[0] << "\n"
                                        << "min: " << hv[1] << "\n"
                                        << "max: " << hv[2] << "\n"
                                        << "isum: " << isum << "\n";
    }

    {
        BL_PROFILE("FabReduce-sum");
        Real hv = 0.0;
        Gpu::DeviceScalar<Real> dv(hv);
        Real* dp = dv.dataPtr();

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            auto const& fab = mf.array(mfi);
            amrex::FabReduce
                (bx, // Box
                 hv, // initial values
            // First lambda works on each cell
            [=] AMREX_GPU_DEVICE (int i, int j, int k, Real* r)  noexcept
            {
                *r += fab(i,j,k);
            },
            // Second lambda does reduce and stores the result in global memory
            // Reduce must be done with Gpu::Reduce[Sum|Min|Max].
            [=] AMREX_GPU_DEVICE (Real r) noexcept
            {
                Gpu::ReduceSum(dp, r);
            });
        }

        hv = dv.dataValue();
        // MPI reduce
        ParallelDescriptor::ReduceRealSum(hv);
        amrex::Print().SetPrecision(17) << "sum: " << hv << "\n";
    }

    {
        BL_PROFILE("FabReduce-min");
        Real hv = std::numeric_limits<Real>::max();
        Gpu::DeviceScalar<Real> dv(hv);
        Real* dp = dv.dataPtr();

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            auto const& fab = mf.array(mfi);
            amrex::FabReduce
                (bx, // Box
                 hv, // initial values
            // First lambda works on each cell
            [=] AMREX_GPU_DEVICE (int i, int j, int k, Real* r)  noexcept
            {
                *r = amrex::min(fab(i,j,k), *r);
            },
            // Second lambda does reduce and stores the result in global memory
            // Reduce must be done with Gpu::Reduce[Sum|Min|Max].
            [=] AMREX_GPU_DEVICE (Real r) noexcept
            {
                Gpu::ReduceMin(dp, r);
            });
        }

        hv = dv.dataValue();
        // MPI reduce
        ParallelDescriptor::ReduceRealMin(hv);
        amrex::Print().SetPrecision(17) << "min: " << hv << "\n";
    }

    {
        BL_PROFILE("FabReduce-max");
        Real hv = std::numeric_limits<Real>::lowest();
        Gpu::DeviceScalar<Real> dv(hv);
        Real* dp = dv.dataPtr();

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            auto const& fab = mf.array(mfi);
            amrex::FabReduce
                (bx, // Box
                 hv, // initial values
            // First lambda works on each cell
            [=] AMREX_GPU_DEVICE (int i, int j, int k, Real* r)  noexcept
            {
                *r = amrex::max(fab(i,j,k), *r);
            },
            // Second lambda does reduce and stores the result in global memory
            // Reduce must be done with Gpu::Reduce[Sum|Min|Max].
            [=] AMREX_GPU_DEVICE (Real r) noexcept
            {
                Gpu::ReduceMax(dp, r);
            });
        }

        hv = dv.dataValue();
        // MPI reduce
        ParallelDescriptor::ReduceRealMax(hv);
        amrex::Print().SetPrecision(17) << "max: " << hv << "\n";
    }

    {
        BL_PROFILE("FabReduce-isum");
        int hv = 0;
        Gpu::DeviceScalar<int> dv(hv);
        int* dp = dv.dataPtr();

        for (MFIter mfi(imf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            auto const& ifab = imf.array(mfi);
            amrex::FabReduce
                (bx, // Box
                 hv, // initial values
            // First lambda works on each cell
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int* r)  noexcept
            {
                *r += ifab(i,j,k);
            },
            // Second lambda does reduce and stores the result in global memory
            // Reduce must be done with Gpu::Reduce[Sum|Min|Max].
            [=] AMREX_GPU_DEVICE (int r) noexcept
            {
                Gpu::ReduceSum(dp, r);
            });
        }

        hv = dv.dataValue();
        // MPI reduce
        ParallelDescriptor::ReduceIntSum(hv);
        amrex::Print() << "isum: " << hv << "\n";
    }

    int N = 1000000;
    Gpu::DeviceVector<Real> vec(N);
    Real* pvec = vec.dataPtr();
    amrex::ParallelFor(N, [=] AMREX_GPU_DEVICE (int i) noexcept {
            pvec[i] = amrex::Random() - 0.5;
    });

    {
        BL_PROFILE("VecReduce");
        Gpu::DeviceScalar<Real> ds(0.0);
        Real* dp = ds.dataPtr();
        amrex::VecReduce(N, // int
                         0.0, // initial value
        // First lambda works on each element
        [=] AMREX_GPU_DEVICE (int i, Real* r) noexcept
        {
            *r += std::abs(pvec[i]);
        },
        // Second lambda does reduce and stores the result in global memory
        // Reduce must be done with Gpu::Reduce[Sum|Min|Max].
        [=] AMREX_GPU_DEVICE (Real r) noexcept
        {
            Gpu::ReduceSum(dp, r);
        });
        // We could do MPI reduce if it is needed.
        amrex::Print().SetPrecision(17) << "1-norm: " << ds.dataValue() << "\n";
    }

    {
        BL_PROFILE("Reduce::Sum");
        Real r = Reduce::Sum(N, pvec);
        amrex::Print().SetPrecision(17) << "Reduce::Sum " << r << "\n";
    }

#ifdef AMREX_USE_GPU
    {
        BL_PROFILE("ThrustReduceSum");
        Real r = thrust::reduce(vec.begin(), vec.end(), 0.0, thrust::plus<Real>());
        amrex::Print().SetPrecision(17) << "thrust::reduce sum " << r << "\n";
    }
#endif


}
