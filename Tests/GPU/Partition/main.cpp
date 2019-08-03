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

    Vector<int> x_host(size);
    Gpu::dtoh_memcpy(x_host.dataPtr(), x.dataPtr(), sizeof(int)*size);

    Gpu::DeviceVector<int> x_amrex(size);
    Gpu::dtod_memcpy(x_amrex.dataPtr(), x.dataPtr(), sizeof(int)*size);

    Gpu::DeviceVector<int> x_amrex_stable(size);
    Gpu::dtod_memcpy(x_amrex_stable.dataPtr(), x.dataPtr(), sizeof(int)*size);

    Gpu::DeviceVector<int> x_thrust(size);
    Gpu::dtod_memcpy(x_thrust.dataPtr(), x.dataPtr(), sizeof(int)*size);

    Vector<int> hx(size);
    Gpu::dtoh_memcpy(hx.dataPtr(), x.dataPtr(), sizeof(int)*size);

    Gpu::synchronize();

    {
        BL_PROFILE("CurrentPartition");
        CurrentPartition(x, [=] AMREX_GPU_DEVICE (int i) -> int {return i % 2 == 0;});
        Gpu::synchronize();
    }

    int neven2;
    {
        BL_PROFILE("amrex::Partition");
        neven2 = amrex::Partition(x_amrex, [=] AMREX_GPU_DEVICE (int i) -> int {return i % 2 == 0;});
        Gpu::synchronize();
    }

    {
        BL_PROFILE("amrex::StablePartition");
        amrex::StablePartition(x_amrex_stable, [=] AMREX_GPU_DEVICE (int i) -> int {return i % 2 == 0;});
        Gpu::synchronize();
    }

    {
        BL_PROFILE("thrust::Partition");
        ThrustPartition(x_thrust, [=] AMREX_GPU_DEVICE (int i) -> int {return i % 2 == 0;});
        Gpu::synchronize();
    }

    {
        // verification
        Vector<int> h(size);
        Vector<int> h_amrex(size);
        Vector<int> h_amrex_stable(size);
        Vector<int> h_thrust(size);
        Gpu::dtoh_memcpy(h.dataPtr(), x.dataPtr(), sizeof(int)*size);
        Gpu::dtoh_memcpy(h_amrex.dataPtr(), x_amrex.dataPtr(), sizeof(int)*size);
        Gpu::dtoh_memcpy(h_amrex_stable.dataPtr(), x_amrex_stable.dataPtr(), sizeof(int)*size);
        Gpu::dtoh_memcpy(h_thrust.dataPtr(), x_thrust.dataPtr(), sizeof(int)*size);

        bool prev = (h[0] % 2 == 0);
        int numevens = prev;
        for (int i = 1; i < size; ++i) {
            bool current = (h[i] % 2 == 0);
            numevens += current;
            if (current != prev && current == true) {
                amrex::Abort("CurrentPartition failed");
            }
            prev = current;
        }

        if (numevens != neven2) {
            amrex::Print() << "CurrentPartition: # of evens = " << numevens << "\n"
                           << "amrex::Partition: # of evens = " << neven2 << "\n";
            amrex::Abort("CurrentPartition or amrex::Partition failed");
        }

        prev = (h_amrex[0] % 2 == 0);
        numevens = prev;
        for (int i = 1; i < size; ++i) {
            bool current = (h_amrex[i] % 2 == 0);
            numevens += current;
            if (current != prev && current == true) {
                amrex::Abort("amrex::Partition failed");
            }
            prev = current;
        }

        if (numevens != neven2) {
            amrex::Abort("amrex::Partition failed");
        }

        Vector<int> x_host_stable = x_host;
        std::stable_partition(x_host_stable.begin(), x_host_stable.end(),
                              [] (int i) -> int { return i % 2 == 0; });
        for (int i = 0; i < size; ++i) {
            if (x_host_stable[i] != h_amrex_stable[i]) {
                amrex::Abort("amrex::StablePartition failed");
            }
        }

#if 0
        amrex::Print() << "--------------------------\n";
        for (int i = 0; i < size; ++i) {
            amrex::Print() << "xxxxx " << hx[i] << ", " << h[i] << ", " << h_amrex[i] << ", " << h_thrust[i] << "\n";
        }
        amrex::Print() << "--------------------------\n";
#endif

#if 0
        prev = (h_thrust[0] % 2 == 0);
        for (int i = 1; i < size; ++i) {
            bool current = (h_thrust[i] % 2 == 0);
            if (current != prev && current == true) {
                amrex::Abort("thrust::Partition failed");
            }
            prev = current;
        }
#endif
    }
}
