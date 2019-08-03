#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_Utility.H>
#include <AMReX_CudaContainers.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Partition.H>

#include <thrust/tuple.h>
#include <thrust/gather.h>
#include <thrust/iterator/zip_iterator.h>

using namespace amrex;

struct tupleAdd
{
    template <typename T>
    AMREX_GPU_HOST_DEVICE
    T operator()(T v1, T v2)  {
        thrust::get<0>(v1) = thrust::get<0>(v2) + thrust::get<0>(v1);
        thrust::get<1>(v1) = thrust::get<1>(v2) + thrust::get<1>(v1);
        return v1 ;
    }
};

void TestPartition();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    TestPartition();
    {
        BL_PROFILE_REGION("TEST");
        TestPartition();
    }
    amrex::Finalize();
}

template <typename T, typename F>
int ThrustPartition (Gpu::DeviceVector<T>& x, F f)
{
    std::size_t N = x.size();
    Gpu::DeviceVector<int> index(N);
    Gpu::DeviceVector<int> output(N);
    Gpu::DeviceVector<T>   x_tmp(N);
    Gpu::DeviceVector<int> func(N);

    thrust::sequence(thrust::device, index.begin(), index.end());

    thrust::transform(thrust::device, x.begin(), x.end(), func.begin(), f);

    auto mid = thrust::partition(thrust::cuda::par(Cuda::The_ThrustCachedAllocator()),
                                 index.begin(), index.end(), func.begin(), f);

    return thrust::distance(index.begin(), mid);
}

template <typename T, typename F>
int CurrentPartition (Gpu::DeviceVector<T>& x, F f)
{
    std::size_t N = x.size();
    Gpu::DeviceVector<int> lo(N);
    Gpu::DeviceVector<int> hi(N);
    Gpu::DeviceVector<int> func(N);
    Gpu::DeviceVector<int> index(N);
    Gpu::DeviceVector<int> output(N);
    Gpu::DeviceVector<T>   x_tmp(N);

    auto index_ptr  = index.dataPtr();
    auto lo_ptr     = lo.dataPtr();
    auto hi_ptr     = hi.dataPtr();
    auto output_ptr = output.dataPtr();
    auto func_ptr   = func.dataPtr();
    auto x_ptr      = x.dataPtr();

    AMREX_FOR_1D ( N, i,
    {
        index_ptr[i] = i;
    
        func_ptr[i] = f(x_ptr[i]);
    
        if (func_ptr[i])
        {
            lo_ptr[i] = 1; hi_ptr[i] = 0;
        }
        else 
        {
            lo_ptr[i] = 0; hi_ptr[i] = 1;
        }
    });

    thrust::exclusive_scan(thrust::cuda::par(Cuda::The_ThrustCachedAllocator()), 
                           thrust::make_zip_iterator(thrust::make_tuple(lo.begin(), hi.begin())),
                           thrust::make_zip_iterator(thrust::make_tuple(lo.end(),   hi.end())),
                           thrust::make_zip_iterator(thrust::make_tuple(lo.begin(), hi.begin())),
                           thrust::tuple<int, int>(0, 0), tupleAdd());
    
    int last_func_host, last_lo_host;
    Gpu::dtoh_memcpy(&last_func_host, func.dataPtr() + N - 1, sizeof(int));
    Gpu::dtoh_memcpy(&last_lo_host, lo.dataPtr() + N - 1, sizeof(int));
    int first_hi = last_func_host + last_lo_host;

    AMREX_FOR_1D ( N, i,
    {
        if (func_ptr[i])
        {
            output_ptr[lo_ptr[i]] = index_ptr[i];
        }
        else 
        {
            output_ptr[hi_ptr[i] + first_hi] = index_ptr[i];
        }
    });
    
    thrust::gather(thrust::device,
                   output.begin(), output.end(),
                   x.dataPtr(), x_tmp.dataPtr());

    x.swap(x_tmp);

    return first_hi;
}

void TestPartition ()
{
    ParmParse pp;
    int size;
    pp.get("size", size);

    constexpr int lower = -1;
    constexpr int upper = 15;

    Gpu::DeviceVector<int> x(size);
    {
        BL_PROFILE("Generate");

        auto x_ptr = x.dataPtr();
        AMREX_PARALLEL_FOR_1D (size, idx,
        {
            x_ptr[idx] = std::ceil(amrex::Random() * (upper - lower + 1)) + lower - 1;
        });
    }

    Gpu::DeviceVector<int> x2(size);
    Gpu::dtod_memcpy(x2.dataPtr(), x.dataPtr(), sizeof(int)*size);

    Gpu::DeviceVector<int> x3(size);
    Gpu::dtod_memcpy(x3.dataPtr(), x.dataPtr(), sizeof(int)*size);

    Gpu::synchronize();

    {
        BL_PROFILE("CurrentPartition");
        CurrentPartition(x, [=] AMREX_GPU_DEVICE (int i) {return i % 2 == 0;});
        Gpu::synchronize();
    }

    {
        BL_PROFILE("amrex::Partition");
        amrex::Partition(x2, [=] AMREX_GPU_DEVICE (int i) {return i % 2 == 0;});
        Gpu::synchronize();
    }

    {
        BL_PROFILE("thrust::Partition");
        ThrustPartition(x3, [=] AMREX_GPU_DEVICE (int i) {return i % 2 == 0;});
        Gpu::synchronize();
    }

    {
        // verification
        Vector<int> h(size);
        Vector<int> h2(size);
        Vector<int> h3(size);
        Gpu::dtoh_memcpy(h.dataPtr(), x.dataPtr(), sizeof(int)*size);
        Gpu::dtoh_memcpy(h2.dataPtr(), x2.dataPtr(), sizeof(int)*size);
        Gpu::dtoh_memcpy(h3.dataPtr(), x3.dataPtr(), sizeof(int)*size);

        bool prev = (h[0] % 2 == 0);
        for (int i = 1; i < size; ++i) {
            bool current = (h[i] % 2 == 0);
            if (current != prev && current == true) {
                amrex::Abort("CurrentPartition failed");
            }
            prev = current;
        }

        prev = (h2[0] % 2 == 0);
        for (int i = 1; i < size; ++i) {
            bool current = (h2[i] % 2 == 0);
            if (current != prev && current == true) {
                amrex::Abort("amrex::Partition failed");
            }
            prev = current;
        }

#if 0
        prev = (h3[0] % 2 == 0);
        for (int i = 1; i < size; ++i) {
            bool current = (h3[i] % 2 == 0);
            if (current != prev && current == true) {
                amrex::Abort("thrust::Partition failed");
            }
            prev = current;
        }
#endif
    }
}
