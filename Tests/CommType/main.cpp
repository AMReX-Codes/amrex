#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuComplex.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    int ret_code = EXIT_SUCCESS;

    {
        int ncells = 64;
        BoxArray ba(Box(IntVect(0), IntVect(ncells-1)));
        ba.maxSize(16);
        ba.convert(IntVect(1));
        DistributionMapping dm(ba);

        constexpr int ncomp = 2;
        IntVect nghost(2);
        Periodicity period{IntVect(ncells)};

        auto value = [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) -> Real
        {
            if (i < 0) {
                i += ncells;
            } else if (i >= ncells) {
                i -= ncells;
            }
            if (j < 0) {
                j += ncells;
            } else if (j >= ncells) {
                j -= ncells;
            }
            if (k < 0) {
                k += ncells;
            } else if (k >= ncells) {
                k -= ncells;
            }
            return n + i*ncomp + j*ncomp*ncells + k*ncomp*ncells*ncells;
        };

        // Test GpuArray
        {
            using T = GpuArray<Real,ncomp>;
            FabArray<BaseFab<T>> fa(ba,dm,1,nghost);
            FabArray<BaseFab<T>> fa2(ba,dm,1,nghost);
            FabArray<BaseFab<T>> fa3(ba,dm,1,nghost);
            auto const& ma = fa.arrays();
            auto const& ma2 = fa2.arrays();
            auto const& ma3 = fa3.arrays();

            ParallelFor(fa, IntVect(0),
                        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                auto const& a = ma[b];
                for (int n = 0; n < ncomp; ++n) {
                    a(i,j,k)[n] = value(i,j,k,n);
                }
            });

            fa.FillBoundary(period);

            fa2.ParallelCopy(fa, 0, 0, 1, IntVect(0), nghost, period);

            fa3.setVal(T{});
            fa3.ParallelAdd(fa, 0, 0, 1, nghost, nghost, period);

            auto mask = OverlapMask(fa3,nghost,period);
            auto const& mma = mask.const_arrays();

            auto err = ParReduce(TypeList<ReduceOpMax,ReduceOpMax,ReduceOpMax>{},
                                 TypeList<Real,Real,Real>{},
                                 fa, nghost,
            [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
                                 -> GpuTuple<Real,Real,Real>
            {
                Real r1 = 0, r2 = 0, r3 = 0;
                auto const& a1 = ma[b];
                auto const& a2 = ma2[b];
                auto const& a3 = ma3[b];
                auto const& m = mma[b];
                for (int n = 0; n < ncomp; ++n) {
                    auto v = value(i,j,k,n);
                    r1 = std::max(r1, std::abs(a1(i,j,k)[n] - v));
                    r2 = std::max(r2, std::abs(a2(i,j,k)[n] - v));
                    r3 = std::max(r3, std::abs(a3(i,j,k)[n] - v*m(i,j,k)));
                }
                return {r1, r2, r3};
            });

            AMREX_ALWAYS_ASSERT(amrex::get<0>(err) == 0);
            AMREX_ALWAYS_ASSERT(amrex::get<1>(err) == 0);
            AMREX_ALWAYS_ASSERT(amrex::get<2>(err) == 0);

            Real errmax = std::max({amrex::get<0>(err),
                                    amrex::get<1>(err),
                                    amrex::get<2>(err)});
            ParallelDescriptor::ReduceRealSum(errmax);
            if (errmax != 0) {
                ret_code = EXIT_FAILURE;
            }
        }

        // Test GpuComplex
        {
            using T = GpuComplex<Real>;
            FabArray<BaseFab<T>> fa(ba,dm,1,nghost);
            FabArray<BaseFab<T>> fa2(ba,dm,1,nghost);
            FabArray<BaseFab<T>> fa3(ba,dm,1,nghost);
            auto const& ma = fa.arrays();
            auto const& ma2 = fa2.arrays();
            auto const& ma3 = fa3.arrays();

            ParallelFor(fa, IntVect(0),
                        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                auto const& a = ma[b];
                a(i,j,k) = T{value(i,j,k,0),value(i,j,k,1)};
            });

            fa.FillBoundary(period);

            fa2.ParallelCopy(fa, 0, 0, 1, IntVect(0), nghost, period);

            fa3.setVal(T{});
            fa3.ParallelAdd(fa, 0, 0, 1, nghost, nghost, period);

            auto mask = OverlapMask(fa3,nghost,period);
            auto const& mma = mask.const_arrays();

            auto err = ParReduce(TypeList<ReduceOpMax,ReduceOpMax,ReduceOpMax>{},
                                 TypeList<Real,Real,Real>{},
                                 fa, nghost,
            [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
                                 -> GpuTuple<Real,Real,Real>
            {
                Real r1 = 0, r2 = 0, r3 = 0;
                auto const& a1 = ma[b];
                auto const& a2 = ma2[b];
                auto const& a3 = ma3[b];
                auto const& m = mma[b];
                auto v = GpuComplex{value(i,j,k,0), value(i,j,k,1)};
                r1 = std::max(r1, amrex::norm(a1(i,j,k) - v));
                r2 = std::max(r2, amrex::norm(a2(i,j,k) - v));
                r3 = std::max(r3, amrex::norm(a3(i,j,k) - v*Real(m(i,j,k))));
                return {r1, r2, r3};
            });

            AMREX_ALWAYS_ASSERT(amrex::get<0>(err) == 0);
            AMREX_ALWAYS_ASSERT(amrex::get<1>(err) == 0);
            AMREX_ALWAYS_ASSERT(amrex::get<2>(err) == 0);

            Real errmax = std::max({amrex::get<0>(err),
                                    amrex::get<1>(err),
                                    amrex::get<2>(err)});
            ParallelDescriptor::ReduceRealSum(errmax);
            if (errmax != 0) {
                ret_code = EXIT_FAILURE;
            }
        }
    }
    amrex::Finalize();

    return ret_code;
}
