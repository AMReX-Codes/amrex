#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>
#include <AMReX_GpuPrint.H>

using namespace amrex;

template <typename T>
struct IsNonnegative {
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    bool operator () (T const val) const noexcept {
        return val >= T(0.);
    }
};

template <typename T>
void test ()
{
    int N = 8000;

    amrex::Gpu::DeviceVector<T> v_d(N);
    auto v_d_ptr = v_d.data();
    amrex::ParallelForRNG(N,
                          [=] AMREX_GPU_DEVICE (int i, RandomEngine const& engine) noexcept
                          {
                              v_d_ptr[i] = T(i);
                          });

    for (int j = 0; j < 1000; ++j) {
        amrex::ParallelFor(N,
                           [=] AMREX_GPU_DEVICE (int i) noexcept
                           {
                               // add -1.0 to v_d_ptr[i] if the result would be non-negative
                               Gpu::Atomic::If(&v_d_ptr[i], T(-1.0), amrex::Plus<T>(), IsNonnegative<T>());
                               // can also use +1.0 amrex::Minus<T>
                               Gpu::Atomic::If(&v_d_ptr[N-i-1], T(1.0), amrex::Minus<T>(), IsNonnegative<T>());
                           });

        amrex::ParallelFor(N,
                           [=] AMREX_GPU_DEVICE (int i) noexcept
                           {
                               // lambdas also work
                               Gpu::Atomic::If(&v_d_ptr[N-i-1], T(1.0), amrex::Minus<T>(),
                                               [=] AMREX_GPU_DEVICE (T tmp) noexcept
                                               {
                                                   return (tmp >= T(0.0));
                                               });
                               Gpu::Atomic::If(&v_d_ptr[i], T(1.0), amrex::Minus<T>(),
                                               [=] AMREX_GPU_DEVICE (T tmp) noexcept
                                               {
                                                   return (tmp >= T(0.0));
                                               });
                           });
    }
    Gpu::streamSynchronize();

    std::vector<T> v_h(N);
    Gpu::copyAsync(Gpu::deviceToHost, v_d.begin(), v_d.end(), v_h.begin());
    Gpu::streamSynchronize();

    // The first 4000 entries should all be 0.0
    for (int i = 0; i < 4000; ++i) {
        AMREX_ALWAYS_ASSERT(v_h[i] == 0);
    }

    // The next 4000 should go from 0 to 3999
    for (int i = 4000, j = 0; i < 8000; ++i, ++j) {
        AMREX_ALWAYS_ASSERT(v_h[i] == j);
    }
}


int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    test<int>();
    test<long>();
    test<float>();
    test<double>();
    test<amrex::Long>();
    amrex::Finalize();
}
