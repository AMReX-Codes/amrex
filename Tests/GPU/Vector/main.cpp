#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>
#include <AMReX_GpuPrint.H>

using namespace amrex;

// v3 should contain the integers from -5 to 5, inclusive
template <template <typename> class Container>
typename std::enable_if<RunOnGpu<typename Container<int>::allocator_type>::value>::type
checkV3 (const Container<int>& c)
{
    const auto c_ptr = c.dataPtr();
    amrex::ParallelFor(11, [=] AMREX_GPU_DEVICE (const int index) noexcept {
            AMREX_ALWAYS_ASSERT(c_ptr[index] == index-5);
        });
    Gpu::Device::streamSynchronize();
}

// v3 should contain the integers from -5 to 5, inclusive
template <template <typename> class Container>
typename std::enable_if<!RunOnGpu<typename Container<int>::allocator_type>::value>::type
checkV3 (const Container<int>& c)
{
    for (int i=-5, index=0; i <= 5; ++i, ++index)
    {
        AMREX_ALWAYS_ASSERT(c[index] == i);
    }
}

template <template <typename> class Container>
void test_container()
{
    Container<int> v = {1, 2, 4};
    v.insert(v.begin(), 0);
    v.insert(v.begin() + 3, 3);
    v.insert(v.end(), 5);
    v.insert(v.begin(), {-2, -1});

    Container<int> v2;
    v2.push_back(-5);
    v2.push_back(-4);
    v2.push_back(-3);
    v.insert(v.begin(), v2.begin(), v2.end());

    v.insert(v.begin(),1,-6);
    v.insert(v.end(), 2, 6);
    v.erase(v.end()-1, v.end());

    v.pop_back();
    v.erase(v.begin(), v.begin()+1);

    Container<int> v3;
    v3.assign(v.begin(), v.end());

    checkV3<Container>(v3);
}

void async_test()
{
    {
        int N = 32*10;  // 10 warps. Arbitrarily chosen.

        Gpu::AsyncVector<Real> vec(N, 0.0);
        auto ptr = vec.dataPtr();

        // Compute-bound test. Should take some time.
        amrex::ParallelFor(N, [=] AMREX_GPU_DEVICE (int n) noexcept
        {
            Real y = ptr[n];
            Real x = 1.0;
            for (int n = 0; n < 20; ++n) {
                Real dx = -(x*x-y) / (2.*x);
                x += dx;
            }
            ptr[n] = x;

            if (n == 0)
#ifdef AMREX_USE_GPU
                { AMREX_DEVICE_PRINTF(" Answer = %1.16f -- should print second.\n", ptr[n]); }
#else
                { std::cout << "Answer = " << ptr[n] << " -- should print first." << std::endl; }
#endif

        });
    }

    // Async Vector now out of scope. Still completes correctly.

#ifdef AMREX_USE_GPU
    amrex::Print() << "Async Synching -- should print first." << std::endl;
#else
    amrex::Print() << "Async Synching -- should print second." << std::endl;
#endif

    Gpu::Device::streamSynchronize();
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        test_container<Gpu::DeviceVector >();
        test_container<Gpu::HostVector   >();
        test_container<Gpu::ManagedVector>();
        test_container<Gpu::PinnedVector> ();
        test_container<Gpu::AsyncVector>  ();

        async_test();

        amrex::Print() << "Passed! \n";
    }
    amrex::Finalize();
}
