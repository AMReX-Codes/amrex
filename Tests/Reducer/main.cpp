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

        // eval(Box, f(iop, i,j,k))
        {
            Reducer<Real> reducer(Vector<ReduceOpType>{
                    ReduceOpType::min, ReduceOpType::max, ReduceOpType::sum,
                    ReduceOpType::sum, ReduceOpType::max});
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                Box const& box = mfi.validbox();
                auto const& a = mf.array(mfi);
                reducer.eval(box, [=] AMREX_GPU_DEVICE (int iop, int i, int j, int k)
                { // 0 <= iop < 5
                    if (iop >= 0 && iop <= 2) { // min, max & sum
                        return a(i,j,k);
                    } else { // 1-norm & inf-norm
                        return std::abs(a(i,j,k));
                    }
                });
            }
            Vector<Real> result = reducer.getResults();

            AMREX_ALWAYS_ASSERT(benchmark1[0] == result[0] &&
                                benchmark1[1] == result[1] &&
                                almostEqual(benchmark1[2], result[2], 100) &&
                                almostEqual(benchmark1[3], result[3], 100) &&
                                benchmark1[4] == result[4]);
        }

        // eval(Box, f(iop, IntVect))
        {
            Reducer<Real> reducer(Vector<ReduceOpType>{
                    ReduceOpType::min, ReduceOpType::max, ReduceOpType::sum,
                    ReduceOpType::sum, ReduceOpType::max});
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                Box const& box = mfi.validbox();
                auto const& a = mf.array(mfi);
                reducer.eval(box, [=] AMREX_GPU_DEVICE (int iop, IntVect const& iv)
                { // 0 <= iop < 5
                    if (iop >= 0 && iop <= 2) { // min, max & sum
                        return a(iv);
                    } else { // 1-norm & inf-norm
                        return std::abs(a(iv));
                    }
                });
            }
            Vector<Real> result = reducer.getResults();

            AMREX_ALWAYS_ASSERT(benchmark1[0] == result[0] &&
                                benchmark1[1] == result[1] &&
                                almostEqual(benchmark1[2], result[2], 100) &&
                                almostEqual(benchmark1[3], result[3], 100) &&
                                benchmark1[4] == result[4]);
        }

        // eval(Box, f(iop, i,j,k,n))
        {
            Reducer<Real> reducer(Vector<ReduceOpType>{
                    ReduceOpType::min, ReduceOpType::max, ReduceOpType::sum,
                    ReduceOpType::sum, ReduceOpType::max});
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                Box const& box = mfi.validbox();
                auto const& a = mf.array(mfi);
                reducer.eval(box, mf.nComp(),
                [=] AMREX_GPU_DEVICE (int iop, int i, int j, int k, int n)
                { // 0 <= iop < 5
                    if (iop >= 0 && iop <= 2) { // min, max & sum
                        return a(i,j,k,n);
                    } else { // 1-norm & inf-norm
                        return std::abs(a(i,j,k,n));
                    }
                });
            }
            Vector<Real> result = reducer.getResults();

            AMREX_ALWAYS_ASSERT(benchmark2[0] == result[0] &&
                                benchmark2[1] == result[1] &&
                                almostEqual(benchmark2[2], result[2], 100) &&
                                almostEqual(benchmark2[3], result[3], 100) &&
                                benchmark2[4] == result[4]);
        }

        // eval(Box, f(iop, IntVect, n))
        {
            Reducer<Real> reducer(Vector<ReduceOpType>{
                    ReduceOpType::min, ReduceOpType::max, ReduceOpType::sum,
                    ReduceOpType::sum, ReduceOpType::max});
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                Box const& box = mfi.validbox();
                auto const& a = mf.array(mfi);
                reducer.eval(box, mf.nComp(),
                [=] AMREX_GPU_DEVICE (int iop, IntVect const& iv, int n)
                { // 0 <= iop < 5
                    if (iop >= 0 && iop <= 2) { // min, max & sum
                        return a(iv,n);
                    } else { // 1-norm & inf-norm
                        return std::abs(a(iv,n));
                    }
                });
            }
            Vector<Real> result = reducer.getResults();

            AMREX_ALWAYS_ASSERT(benchmark2[0] == result[0] &&
                                benchmark2[1] == result[1] &&
                                almostEqual(benchmark2[2], result[2], 100) &&
                                almostEqual(benchmark2[3], result[3], 100) &&
                                benchmark2[4] == result[4]);
        }

        // eval(MultiFab, f(iop, b, i,j,k))
        {
            Reducer<Real> reducer(Vector<ReduceOpType>{
                    ReduceOpType::min, ReduceOpType::max, ReduceOpType::sum,
                    ReduceOpType::sum, ReduceOpType::max});
            auto const& ma = mf.const_arrays();
            reducer.eval(mf, IntVect(0),
                         [=] AMREX_GPU_DEVICE (int iop, int b, int i, int j, int k)
            { // 0 <= iop < 5
                if (iop >= 0 && iop <= 2) { // min, max & sum
                    return ma[b](i,j,k);
                } else { // 1-norm & inf-norm
                    return std::abs(ma[b](i,j,k));
                }
            });
            Vector<Real> result = reducer.getResults();

            AMREX_ALWAYS_ASSERT(benchmark1[0] == result[0] &&
                                benchmark1[1] == result[1] &&
                                almostEqual(benchmark1[2], result[2], 100) &&
                                almostEqual(benchmark1[3], result[3], 100) &&
                                benchmark1[4] == result[4]);
        }

        // eval(MultiFab, f(iop, b, i,j,k,n))
        {
            Reducer<Real> reducer(Vector<ReduceOpType>{
                    ReduceOpType::min, ReduceOpType::max, ReduceOpType::sum,
                    ReduceOpType::sum, ReduceOpType::max});
            auto const& ma = mf.const_arrays();
            reducer.eval(mf, IntVect(0), mf.nComp(),
                [=] AMREX_GPU_DEVICE (int iop, int b, int i, int j, int k, int n)
            { // 0 <= iop < 5
                if (iop >= 0 && iop <= 2) { // min, max & sum
                    return ma[b](i,j,k,n);
                } else { // 1-norm & inf-norm
                    return std::abs(ma[b](i,j,k,n));
                }
            });
            Vector<Real> result = reducer.getResults();

            AMREX_ALWAYS_ASSERT(benchmark2[0] == result[0] &&
                                benchmark2[1] == result[1] &&
                                almostEqual(benchmark2[2], result[2], 100) &&
                                almostEqual(benchmark2[3], result[3], 100) &&
                                benchmark2[4] == result[4]);
        }

        // eval(int, f(iop, i))
        {
            Reducer<Real> reducer(Vector<ReduceOpType>{
                    ReduceOpType::min, ReduceOpType::max, ReduceOpType::sum,
                    ReduceOpType::sum, ReduceOpType::max});
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                auto const* p = mf[mfi].dataPtr();
                Long np = mfi.fabbox().numPts();
                reducer.eval(np, [=] AMREX_GPU_DEVICE (int iop, Long i)
                { // 0 <= iop < 5
                    if (iop >= 0 && iop <= 2) { // min, max & sum
                        return p[i];
                    } else { // 1-norm & inf-norm
                        return std::abs(p[i]);
                    }
                });
            }
            Vector<Real> result = reducer.getResults();

            AMREX_ALWAYS_ASSERT(benchmark1[0] == result[0] &&
                                benchmark1[1] == result[1] &&
                                almostEqual(benchmark1[2], result[2], 100) &&
                                almostEqual(benchmark1[3], result[3], 100) &&
                                benchmark1[4] == result[4]);
        }
    }

    amrex::Finalize();
}
