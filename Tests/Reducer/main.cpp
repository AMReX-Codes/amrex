#include <AMReX.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Reduce.H>

using namespace amrex;


int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    int n_cell = 64;
    int max_grid_size = 32;
    {
        ParmParse pp;
        pp.query("n_cell", n_cell);
        pp.query("max_grid_size", max_grid_size);
    }

    {
        BoxArray ba(Box(IntVect(0),IntVect(n_cell-1)));
        ba.maxSize(max_grid_size);
        MultiFab mf(ba, DistributionMapping{ba}, 2, 0);
        FillRandom(mf, 0, 2);
        mf.plus(Real(-0.2), 0, 2);

        // No need to use MPI in testing local reduce. Hence true.
        Vector<Real> benchmark1{mf.min(0,0,true), mf.max(0,0,true),
            mf.sum(0,true), mf.norm1(0,0,true), mf.norminf(0,0,true)};

        Vector<Real> benchmark2{mf.min(1,0,true), mf.max(1,0,true),
            mf.sum(1,true), mf.norm1(1,0,true), mf.norminf(1,0,true)};

        benchmark2[0] = std::min(benchmark2[0], benchmark1[0]);
        benchmark2[1] = std::max(benchmark2[1], benchmark1[1]);
        benchmark2[2] += benchmark1[2];
        benchmark2[3] += benchmark1[3];
        benchmark2[4] = std::max(benchmark2[4], benchmark1[4]);

        // one reduce
        {
            Reducer<ReduceOpMin,Real> reducer;
            auto const& ma = mf.const_arrays();
            reducer.eval(mf, IntVect(0), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                return ma[b](i,j,k);
            });
            auto result = reducer.getResult();
            AMREX_ALWAYS_ASSERT(benchmark1[0] == amrex::get<0>(result));
        }

        // eval(Box, f(i,j,k))
        {
            Reducer<TypeList<ReduceOpMin,ReduceOpMax,ReduceOpSum,ReduceOpSum,ReduceOpMax>,
                    TypeMultiplier<TypeList,Real[5]>> reducer;
            using Result_t = typename decltype(reducer)::Result_t;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                Box const& box = mfi.validbox();
                auto const& a = mf.array(mfi);
                reducer.eval(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> Result_t
                {
                    auto v = a(i,j,k);
                    auto vabs = std::abs(v);
                    return {v,v,v,vabs,vabs};
                });
            }
            auto result = amrex::tupleToArray(reducer.getResult());

            AMREX_ALWAYS_ASSERT(benchmark1[0] == result[0] &&
                                benchmark1[1] == result[1] &&
                                almostEqual(benchmark1[2], result[2], 100) &&
                                almostEqual(benchmark1[3], result[3], 100) &&
                                benchmark1[4] == result[4]);
        }

        // eval(Box, f(IntVect))
        {
            Reducer<TypeList<ReduceOpMin,ReduceOpMax,ReduceOpSum,ReduceOpSum,ReduceOpMax>,
                    TypeList<Real,Real,Real,Real,Real>> reducer;
            using Result_t = typename decltype(reducer)::Result_t;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                Box const& box = mfi.validbox();
                auto const& a = mf.array(mfi);
                reducer.eval(box, [=] AMREX_GPU_DEVICE (IntVect const& iv) -> Result_t
                {
                    auto v = a(iv);
                    auto vabs = std::abs(v);
                    return {v,v,v,vabs,vabs};
                });
            }
            auto result = amrex::tupleToArray(reducer.getResult());

            AMREX_ALWAYS_ASSERT(benchmark1[0] == result[0] &&
                                benchmark1[1] == result[1] &&
                                almostEqual(benchmark1[2], result[2], 100) &&
                                almostEqual(benchmark1[3], result[3], 100) &&
                                benchmark1[4] == result[4]);
        }

        // eval(Box, f(i,j,k,n))
        {
            Reducer<TypeList<ReduceOpMin,ReduceOpMax,ReduceOpSum,ReduceOpSum,ReduceOpMax>,
                    TypeMultiplier<TypeList,Real[5]>> reducer;
            using Result_t = typename decltype(reducer)::Result_t;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                Box const& box = mfi.validbox();
                auto const& a = mf.array(mfi);
                reducer.eval(box, mf.nComp(),
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) -> Result_t
                {
                    auto v = a(i,j,k,n);
                    auto vabs = std::abs(v);
                    return {v,v,v,vabs,vabs};
                });
            }
            auto result = amrex::tupleToArray(reducer.getResult());

            AMREX_ALWAYS_ASSERT(benchmark2[0] == result[0] &&
                                benchmark2[1] == result[1] &&
                                almostEqual(benchmark2[2], result[2], 100) &&
                                almostEqual(benchmark2[3], result[3], 100) &&
                                benchmark2[4] == result[4]);
        }

        // eval(Box, f(IntVect, n))
        {
            Reducer<TypeList<ReduceOpMin,ReduceOpMax,ReduceOpSum,ReduceOpSum,ReduceOpMax>,
                    TypeMultiplier<TypeList,Real[5]>> reducer;
            using Result_t = typename decltype(reducer)::Result_t;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                Box const& box = mfi.validbox();
                auto const& a = mf.array(mfi);
                reducer.eval(box, mf.nComp(),
                [=] AMREX_GPU_DEVICE (IntVect const& iv, int n) -> Result_t
                {
                    auto v = a(iv,n);
                    auto vabs = std::abs(v);
                    return {v,v,v,vabs,vabs};
                });
            }
            auto result = amrex::tupleToArray(reducer.getResult());

            AMREX_ALWAYS_ASSERT(benchmark2[0] == result[0] &&
                                benchmark2[1] == result[1] &&
                                almostEqual(benchmark2[2], result[2], 100) &&
                                almostEqual(benchmark2[3], result[3], 100) &&
                                benchmark2[4] == result[4]);
        }

        // eval(MultiFab, f(b, i,j,k))
        {
            Reducer<TypeList<ReduceOpMin,ReduceOpMax,ReduceOpSum,ReduceOpSum,ReduceOpMax>,
                    TypeMultiplier<TypeList,Real,Real[3],Real>> reducer;
            using Result_t = typename decltype(reducer)::Result_t;
            auto const& ma = mf.const_arrays();
            reducer.eval(mf, IntVect(0),
                         [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) -> Result_t
            {
                auto v = ma[b](i,j,k);
                auto vabs = std::abs(v);
                return {v,v,v,vabs,vabs};
            });
            auto result = amrex::tupleToArray(reducer.getResult());

            AMREX_ALWAYS_ASSERT(benchmark1[0] == result[0] &&
                                benchmark1[1] == result[1] &&
                                almostEqual(benchmark1[2], result[2], 100) &&
                                almostEqual(benchmark1[3], result[3], 100) &&
                                benchmark1[4] == result[4]);
        }

        // eval(MultiFab, f(b, i,j,k,n))
        {
            Reducer<TypeList<ReduceOpMin,ReduceOpMax,ReduceOpSum,ReduceOpSum,ReduceOpMax>,
                    TypeMultiplier<TypeList,Real,Real[3],Real>> reducer;
            using Result_t = typename decltype(reducer)::Result_t;
            auto const& ma = mf.const_arrays();
            reducer.eval(mf, IntVect(0), mf.nComp(),
                [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n) -> Result_t
            {
                auto v = ma[b](i,j,k,n);
                auto vabs = std::abs(v);
                return {v,v,v,vabs,vabs};
            });
            auto result = amrex::tupleToArray(reducer.getResult());

            AMREX_ALWAYS_ASSERT(benchmark2[0] == result[0] &&
                                benchmark2[1] == result[1] &&
                                almostEqual(benchmark2[2], result[2], 100) &&
                                almostEqual(benchmark2[3], result[3], 100) &&
                                benchmark2[4] == result[4]);
        }

        // eval(int, f(i))
        {
            Reducer<TypeList<ReduceOpMin,ReduceOpMax,ReduceOpSum,ReduceOpSum,ReduceOpMax>,
                    TypeMultiplier<TypeList,Real[5]>> reducer;
            using Result_t = typename decltype(reducer)::Result_t;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                auto const* p = mf[mfi].dataPtr();
                Long np = mfi.fabbox().numPts();
                reducer.eval(np, [=] AMREX_GPU_DEVICE (Long i) -> Result_t
                {
                    auto v = p[i];
                    auto vabs = std::abs(v);
                    return {v,v,v,vabs,vabs};
                });
            }
            auto result = amrex::tupleToArray(reducer.getResult());

            AMREX_ALWAYS_ASSERT(benchmark1[0] == result[0] &&
                                benchmark1[1] == result[1] &&
                                almostEqual(benchmark1[2], result[2], 100) &&
                                almostEqual(benchmark1[3], result[3], 100) &&
                                benchmark1[4] == result[4]);
        }
    }

    amrex::Finalize();
}
