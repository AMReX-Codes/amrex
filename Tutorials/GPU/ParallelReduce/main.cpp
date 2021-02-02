
#include <AMReX.H>
#include <AMReX_Utility.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

//#ifdef AMREX_USE_HIP
//#include <thrust/reduce.h>
//#endif

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
    int nghost = 0;
    {
        ParmParse pp;
        pp.query("ncell", ncell);
        pp.query("max_grid_size", max_grid_size);
        pp.query("nghost", nghost);
    }

    BoxArray ba;
    {
        Box domain_box(IntVect(0), IntVect(ncell-1));
        ba.define(domain_box);
        ba.maxSize(max_grid_size);
    }

    MultiFab mf(ba,DistributionMapping{ba},1,nghost);
    iMultiFab imf(ba,mf.DistributionMap(),1,nghost);
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.fabbox();
        auto const& fab = mf.array(mfi);
        auto const& ifab = imf.array(mfi);

        amrex::ParallelForRNG(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, RandomEngine const& engine) noexcept
        {
            fab(i,j,k) = amrex::Random(engine);
            ifab(i,j,k) = amrex::Random_int(2,engine);
        });
    }

    int N = 1000000;
    Gpu::DeviceVector<Real> vec(N);
    Real* pvec = vec.dataPtr();
    amrex::ParallelForRNG( N,
    [=] AMREX_GPU_DEVICE (int i, RandomEngine const& engine) noexcept
    {
        pvec[i] = amrex::Random(engine) - 0.5;
    });

    {
        BL_PROFILE("ParallelForReduction-box-3");

        Gpu::Buffer<Real> da({0.0, std::numeric_limits<Real>::max(),
                                       std::numeric_limits<Real>::lowest()});
        Real* dp = da.data();
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.fabbox();
            Array4<Real const> const& fab = mf.const_array(mfi);
            amrex::ParallelFor(Gpu::KernelInfo().setReduction(true), bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, Gpu::Handler const& handler) noexcept
            {
                Gpu::deviceReduceSum(dp  , fab(i,j,k), handler);
                Gpu::deviceReduceMin(dp+1, fab(i,j,k), handler);
                Gpu::deviceReduceMax(dp+2, fab(i,j,k), handler);
            });
        }
        Real* hp = da.copyToHost();
        ParallelDescriptor::ReduceRealSum(hp[0]);
        ParallelDescriptor::ReduceRealMin(hp[1]);
        ParallelDescriptor::ReduceRealMax(hp[2]);
        amrex::Print().SetPrecision(17) << "sum: "  << hp[0] << "\n"
                                        << "min: "  << hp[1] << "\n"
                                        << "max: "  << hp[2] << "\n";
    }

   {
        BL_PROFILE("ParallelForReduction-box-sum");

        Gpu::Buffer<Real> da({0.0});
        Real* dp = da.data();
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.fabbox();
            Array4<Real const> const& fab = mf.const_array(mfi);
            amrex::ParallelFor(Gpu::KernelInfo().setReduction(true), bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, Gpu::Handler const& handler) noexcept
            {
                Gpu::deviceReduceSum(dp, fab(i,j,k), handler);
            });
        }
        Real* hp = da.copyToHost();
        ParallelDescriptor::ReduceRealSum(hp[0]);
        amrex::Print().SetPrecision(17) << "sum: "  << hp[0] << "\n";
    }

   {
        BL_PROFILE("ParallelForReduction-box-isum");

        Gpu::Buffer<Long> da({0});
        Long* dp = da.data();
        for (MFIter mfi(imf); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.fabbox();
            Array4<int const> const& ifab = imf.const_array(mfi);
            amrex::ParallelFor(Gpu::KernelInfo().setReduction(true), bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, Gpu::Handler const& handler) noexcept
            {
                Gpu::deviceReduceSum<Long>(dp, ifab(i,j,k), handler);
            });
        }
        Long* hp = da.copyToHost();
        ParallelDescriptor::ReduceLongSum(hp[0]);
        amrex::Print().SetPrecision(17) << "isum: "  << hp[0] << "\n";
    }

    {
        BL_PROFILE("ParallelForReduction-vec-1");
        Gpu::Buffer<Real> da({0.0});
        Real* dp = da.data();
        amrex::ParallelFor(Gpu::KernelInfo().setReduction(true), N,
        [=] AMREX_GPU_DEVICE (int i, Gpu::Handler const& handler) noexcept
        {
            Gpu::deviceReduceSum(dp, amrex::Math::abs(pvec[i]), handler);
        });
        Real* hp = da.copyToHost();
        ParallelDescriptor::ReduceRealSum(hp[0]);
        amrex::Print().SetPrecision(17) << "1-norm: "  << hp[0] << "\n";
    }

    {
        BL_PROFILE("MultiFab::sum");
        amrex::Print().SetPrecision(17) << "sum: " << mf.sum() << "\n";
    }

    {
        BL_PROFILE("MultiFab::min");
        amrex::Print().SetPrecision(17) << "min: " << mf.min(0, nghost) << "\n";
    }

    {
        BL_PROFILE("MultiFab::max");
        amrex::Print().SetPrecision(17) << "max: " << mf.max(0, nghost) << "\n";
    }

    {
        BL_PROFILE("iMultiFab::sum");
        amrex::Print() << "isum: " << imf.sum(0) << "\n";
    }

    {
        BL_PROFILE("FabReduceTuple");

        ReduceOps<ReduceOpSum, ReduceOpMin, ReduceOpMax, ReduceOpSum> reduce_op;
        ReduceData<Real, Real, Real, Long> reduce_data(reduce_op);
        // For iMultiFab::sum, we use Long to avoid overflow.
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.fabbox();
            auto const& fab = mf.array(mfi);
            auto const& ifab = imf.array(mfi);
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                Real x =  fab(i,j,k);
                Long ix = static_cast<Long>(ifab(i,j,k));
                return {x,x,x,ix};
            });
        }

        ReduceTuple hv = reduce_data.value();
        // MPI reduce
        ParallelDescriptor::ReduceRealSum(amrex::get<0>(hv));
        ParallelDescriptor::ReduceRealMin(amrex::get<1>(hv));
        ParallelDescriptor::ReduceRealMax(amrex::get<2>(hv));
        ParallelDescriptor::ReduceLongSum(amrex::get<3>(hv));
        amrex::Print().SetPrecision(17) << "sum: "  << get<0>(hv) << "\n"
                                        << "min: "  << get<1>(hv) << "\n"
                                        << "max: "  << get<2>(hv) << "\n"
                                        << "isum: " << get<3>(hv) << "\n";
    }

    {
        BL_PROFILE("FabReduceTuple-box");

        ReduceOps<ReduceOpSum, ReduceOpMin, ReduceOpMax, ReduceOpSum> reduce_op;
        ReduceData<Real, Real, Real, Long> reduce_data(reduce_op);
        // For iMultiFab::sum, we use Long to avoid overflow.
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.fabbox();
            auto const& fab = mf.array(mfi);
            auto const& ifab = imf.array(mfi);
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (Box const& b) -> ReduceTuple
            {
                Real rsum = 0.;
                Real rmin =  1.e30; // If not because of cuda 9.2,
                Real rmax = -1.e30; // we should use numeric_limits.
                Long lsum = 0;
                amrex::Loop(b,
                [=,&rsum,&rmin,&rmax,&lsum] (int i, int j, int k) {
                    Real x =  fab(i,j,k);
                    Long ix = static_cast<Long>(ifab(i,j,k));
                    rsum += x;
                    rmin = amrex::min(rmin,x);
                    rmax = amrex::max(rmax,x);
                    lsum += ix;
                });
                return {rsum,rmin,rmax,lsum};
            });
        }

        ReduceTuple hv = reduce_data.value();
        // MPI reduce
        ParallelDescriptor::ReduceRealSum(amrex::get<0>(hv));
        ParallelDescriptor::ReduceRealMin(amrex::get<1>(hv));
        ParallelDescriptor::ReduceRealMax(amrex::get<2>(hv));
        ParallelDescriptor::ReduceLongSum(amrex::get<3>(hv));
        amrex::Print().SetPrecision(17) << "sum: "  << get<0>(hv) << "\n"
                                        << "min: "  << get<1>(hv) << "\n"
                                        << "max: "  << get<2>(hv) << "\n"
                                        << "isum: " << get<3>(hv) << "\n";
    }

    {
        BL_PROFILE("FabReduce-sum");

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.fabbox();
            auto const& fab = mf.array(mfi);
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                return fab(i,j,k);
            });
        }

        Real hv = amrex::get<0>(reduce_data.value());
        // MPI reduce
        ParallelDescriptor::ReduceRealSum(hv);
        amrex::Print().SetPrecision(17) << "sum: " << hv << "\n";
    }

    {
        BL_PROFILE("FabReduce-min");

        ReduceOps<ReduceOpMin> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.fabbox();
            auto const& fab = mf.array(mfi);
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                return fab(i,j,k);
            });
        }

        Real hv = amrex::get<0>(reduce_data.value());
        // MPI reduce
        ParallelDescriptor::ReduceRealMin(hv);
        amrex::Print().SetPrecision(17) << "min: " << hv << "\n";
    }

    {
        BL_PROFILE("FabReduce-max");

        ReduceOps<ReduceOpMax> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.fabbox();
            auto const& fab = mf.array(mfi);
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                return fab(i,j,k);
            });
        }

        Real hv = amrex::get<0>(reduce_data.value());
        // MPI reduce
        ParallelDescriptor::ReduceRealMax(hv);
        amrex::Print().SetPrecision(17) << "max: " << hv << "\n";
    }

    {
        BL_PROFILE("FabReduce-isum-long");

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Long> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(imf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.fabbox();
            auto const& ifab = imf.array(mfi);
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                return static_cast<Long>(ifab(i,j,k));
            });
        }

        Long hv = amrex::get<0>(reduce_data.value());
        // MPI reduce
        ParallelDescriptor::ReduceLongSum(hv);
        amrex::Print() << "isum: " << hv << "\n";
    }

    {
        // Changing types to take advantage of available hardware acceleration.
        // Recommeded version for GPUs. (~60x faster).
        BL_PROFILE("FabReduce-isum-unsigned-long-long");

        long long points = 0;
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<unsigned long long> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(imf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            auto const& ifab = imf.array(mfi);
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                return { static_cast<unsigned long long>(ifab(i,j,k)) -
                         static_cast<unsigned long long>(INT_MIN) };
            });
            points += bx.numPts();
        }

        Long hv = static_cast<Long>( static_cast<long long>(amrex::get<0>(reduce_data.value()))
                                     + static_cast<long long>(INT_MIN)*points );
        // MPI reduce
        ParallelDescriptor::ReduceLongSum(hv);
        amrex::Print() << "isum: " << hv << "\n";
    }

    {
        BL_PROFILE("VecReduce");

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        reduce_op.eval(N, reduce_data,
        [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
        {
            return amrex::Math::abs(pvec[i]);
        });

        Real hv = amrex::get<0>(reduce_data.value());
        // We could do MPI reduce if it is needed.
        amrex::Print().SetPrecision(17) << "1-norm: " << hv << "\n";
    }

    {
        BL_PROFILE("Reduce::Sum");
        Real r = Reduce::Sum(N, pvec);
        amrex::Print().SetPrecision(17) << "Reduce::Sum " << r << "\n";
    }

    {
        BL_PROFILE("Reduce::Min");
        Real r = Reduce::Min(N, pvec);
        amrex::Print().SetPrecision(17) << "Reduce::Min " << r << "\n";
    }

    {
        BL_PROFILE("Reduce::Max");
        Real r = Reduce::Max(N, pvec);
        amrex::Print().SetPrecision(17) << "Reduce::Max " << r << "\n";
    }

    {
        BL_PROFILE("Reduce::MinMax");
        std::pair<Real,Real> r = Reduce::MinMax(N, pvec);
        amrex::Print().SetPrecision(17) << "Reduce::MinMax " << r.first
                                        << ", " << r.second << "\n";
    }

#if 0
#if defined(AMREX_USE_GPU)
    {
        BL_PROFILE("ThrustReduceSum");
        Real r = thrust::reduce(thrust::device, vec.begin(), vec.end(), 0.0, thrust::plus<Real>());
        amrex::Print().SetPrecision(17) << "thrust::reduce sum " << r << "\n";
    }
#endif
#endif
}
