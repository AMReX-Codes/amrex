
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

        ReduceOps<ReduceOpSum, ReduceOpMin, ReduceOpMax, ReduceOpSum> reduce_op;
        ReduceData<Real, Real, Real, long> reduce_data(reduce_op);
        // For iMultiFab::sum, we use long to avoid overflow.
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            auto const& fab = mf.array(mfi);
            auto const& ifab = imf.array(mfi);
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                Real x =  fab(i,j,k);
                long ix = static_cast<long>(ifab(i,j,k));
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
        ReduceData<Real, Real, Real, long> reduce_data(reduce_op);
        // For iMultiFab::sum, we use long to avoid overflow.
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            auto const& fab = mf.array(mfi);
            auto const& ifab = imf.array(mfi);
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (Box const& bx) -> ReduceTuple
            {
                Real rsum = 0.;
                Real rmin =  1.e30; // If not because of cuda 9.2,
                Real rmax = -1.e30; // we should use numeric_limits.
                long lsum = 0;
                amrex::Loop(bx,
                [=,&rsum,&rmin,&rmax,&lsum] (int i, int j, int k) noexcept {
                    Real x =  fab(i,j,k);
                    long ix = static_cast<long>(ifab(i,j,k));
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
            const Box& bx = mfi.validbox();
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
            const Box& bx = mfi.validbox();
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
            const Box& bx = mfi.validbox();
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
        BL_PROFILE("FabReduce-isum");

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<long> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(imf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            auto const& ifab = imf.array(mfi);
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                return static_cast<long>(ifab(i,j,k));
            });
        }

        long hv = amrex::get<0>(reduce_data.value());
        // MPI reduce
        ParallelDescriptor::ReduceLongSum(hv);
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

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        reduce_op.eval(N, reduce_data,
        [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
        {
            return std::abs(pvec[i]);
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

#ifdef AMREX_USE_GPU
    {
        BL_PROFILE("ThrustReduceSum");
        Real r = thrust::reduce(vec.begin(), vec.end(), 0.0, thrust::plus<Real>());
        amrex::Print().SetPrecision(17) << "thrust::reduce sum " << r << "\n";
    }
#endif


}
